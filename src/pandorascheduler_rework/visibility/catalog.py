"""Public API for generating visibility artifacts."""

from __future__ import annotations

import concurrent.futures
import io
import logging
import multiprocessing
import os
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.time import Time
from tqdm import tqdm

import pyarrow as pa
import pyarrow.parquet as pq

from pandorascheduler_rework.config import PandoraSchedulerConfig, resolve_data_subdir
from pandorascheduler_rework.utils.io import read_csv_cached, read_parquet_cached

from .constraints import (
    _R_EARTH_KM,
    _normalise,
    compute_visibility_with_constraints,
    detect_orbit_boundaries,
    orbit_slices_from_boundaries,
)
from .geometry import (
    MinuteCadence,
    build_minute_cadence,
    compute_saa_crossings,
    interpolate_gmat_ephemeris,
)

LOGGER = logging.getLogger(__name__)


def _write_visibility_parquet(
    df: pd.DataFrame,
    path_or_buf: Path | io.BytesIO,
    config: PandoraSchedulerConfig,
) -> None:
    """Write *df* to parquet with keepout-angle metadata in the schema."""
    table = pa.Table.from_pandas(df, preserve_index=False)
    existing = table.schema.metadata or {}
    existing.update(
        {
            b"pandora.visibility_sun_deg": str(config.sun_avoidance_deg).encode(),
            b"pandora.visibility_moon_deg": str(config.moon_avoidance_deg).encode(),
            b"pandora.visibility_earth_deg": str(config.earth_avoidance_deg).encode(),
        }
    )
    table = table.replace_schema_metadata(existing)
    pq.write_table(
        table,
        path_or_buf,
        compression="snappy",
        write_statistics=False,
        use_dictionary=False,
    )


# ---------------------------------------------------------------------------
# Worker state for multiprocessing (set once per worker via _init_worker)
# ---------------------------------------------------------------------------
_worker_payload: dict | None = None
_worker_config: PandoraSchedulerConfig | None = None


def _init_worker(
    payload: dict,
    config: PandoraSchedulerConfig,
) -> None:
    """Initialise per-worker shared state (called once when worker starts)."""
    global _worker_payload, _worker_config  # noqa: PLW0603
    _worker_payload = payload
    _worker_config = config


def _worker_build_star(
    star_name: str,
    star_coord: SkyCoord,
    is_exoplanet: bool,
) -> tuple[str, bytes]:
    """Build visibility for one star and return parquet bytes.

    Runs inside a worker process.  Reads *_worker_payload* and
    *_worker_config* set by :func:`_init_worker`.
    """
    assert _worker_payload is not None and _worker_config is not None
    df = _build_star_visibility(_worker_payload, star_coord, _worker_config)
    if not is_exoplanet:
        df["Time(MJD_UTC)"] = np.round(df["Time(MJD_UTC)"], 6)
    buf = io.BytesIO()
    _write_visibility_parquet(df, buf, _worker_config)
    return star_name, buf.getvalue()


def build_visibility_catalog(
    config: PandoraSchedulerConfig,
    target_list: Path,
    partner_list: Path | None = None,
    output_subpath: str = "targets",
) -> None:
    """Generate visibility outputs for the requested targets."""

    if not config.output_dir:
        raise ValueError("config.output_dir is required for visibility generation")

    data_subdir = resolve_data_subdir(
        config.extra_inputs,
        sun_avoidance_deg=config.sun_avoidance_deg,
        moon_avoidance_deg=config.moon_avoidance_deg,
        earth_avoidance_deg=config.earth_avoidance_deg,
        earth_avoidance_day_deg=config.earth_avoidance_day_deg,
    )
    output_root = config.output_dir / data_subdir / output_subpath
    output_root.mkdir(parents=True, exist_ok=True)

    target_path = target_list if target_list.is_absolute() else target_list.resolve()
    gmat_path = (
        config.gmat_ephemeris
        if config.gmat_ephemeris.is_absolute()
        else config.gmat_ephemeris.resolve()
    )

    target_manifest = _load_target_manifest(target_path, config.target_filters)
    if target_manifest.empty:
        LOGGER.info("No targets matched visibility configuration; skipping build.")
        return

    # Check if any star visibility files need to be generated
    stars_to_generate = []
    for _, row in target_manifest.iterrows():
        star_name = str(row.get("Star Name", ""))
        output_path = output_root / star_name / f"Visibility for {star_name}.parquet"
        if not output_path.exists() or config.force_regenerate:
            stars_to_generate.append((star_name, row))

    # Only compute expensive ephemeris/payload if we need to generate files
    if stars_to_generate:
        cadence = build_minute_cadence(config.window_start, config.window_end)
        ephemeris = interpolate_gmat_ephemeris(gmat_path, cadence)
        base_payload = _build_base_payload(ephemeris, cadence)
        star_metadata = _build_star_metadata(target_manifest)

        is_exoplanet = "exoplanet" in target_path.name.lower()

        # Resolve coordinates up front (fast, needs star_metadata)
        work_items: list[tuple[str, SkyCoord]] = []
        for star_name, row in stars_to_generate:
            star_coord = _resolve_star_coord(row, star_metadata)
            work_items.append((star_name, star_coord))

        n_stars = len(work_items)
        max_workers = config.parallel_workers or (os.cpu_count() or 1)
        n_workers = min(n_stars, max_workers)

        if n_workers > 1:
            LOGGER.info(
                "Generating visibility for %d stars using %d workers",
                n_stars,
                n_workers,
            )
            available_methods = multiprocessing.get_all_start_methods()
            start_method = (
                "forkserver" if "forkserver" in available_methods else "spawn"
            )
            # Prefer forkserver where available to avoid macOS fork-safety
            # warnings; otherwise fall back to spawn for portability.
            ctx = multiprocessing.get_context(start_method)
            with concurrent.futures.ProcessPoolExecutor(
                max_workers=n_workers,
                mp_context=ctx,
                initializer=_init_worker,
                initargs=(base_payload, config),
            ) as executor:
                futures = {
                    executor.submit(
                        _worker_build_star, name, coord, is_exoplanet
                    ): name
                    for name, coord in work_items
                }
                progress = tqdm(
                    concurrent.futures.as_completed(futures),
                    total=n_stars,
                    desc="Visibility",
                    unit="star",
                    disable=not config.show_progress,
                )
                for future in progress:
                    star_name, parquet_bytes = future.result()
                    output_dir = output_root / star_name
                    output_dir.mkdir(parents=True, exist_ok=True)
                    out_path = (
                        output_dir / f"Visibility for {star_name}.parquet"
                    )
                    out_path.write_bytes(parquet_bytes)
        else:
            LOGGER.info(
                "Generating visibility for %d star(s) (serial)", n_stars
            )
            for star_name, star_coord in tqdm(
                work_items,
                desc="Visibility",
                unit="star",
                disable=not config.show_progress,
            ):
                output_dir = output_root / star_name
                output_dir.mkdir(parents=True, exist_ok=True)
                output_path = (
                    output_dir / f"Visibility for {star_name}.parquet"
                )

                visibility_df = _build_star_visibility(
                    base_payload, star_coord, config
                )

                if not is_exoplanet:
                    visibility_df["Time(MJD_UTC)"] = np.round(
                        visibility_df["Time(MJD_UTC)"], 6
                    )

                _write_visibility_parquet(visibility_df, output_path, config)
    else:
        # Still need star_metadata for planet transits
        star_metadata = _build_star_metadata(target_manifest)

    planet_manifests: list[tuple[pd.DataFrame, Path]] = [(target_manifest, target_path)]

    if partner_list is not None:
        partner_path = (
            partner_list if partner_list.is_absolute() else partner_list.resolve()
        )
        partner_manifest = _load_target_manifest(partner_path, config.target_filters)
        if not partner_manifest.empty:
            star_metadata.update(_build_star_metadata(partner_manifest))
            planet_manifests.append((partner_manifest, partner_path))

    generated_planets: list[tuple[str, str]] = []
    for manifest, manifest_path in planet_manifests:
        generated_planets.extend(
            _build_planet_transits(
                manifest,
                manifest_path,
                output_root,
                star_metadata,
                config,
            )
        )

    _apply_transit_overlaps(generated_planets, output_root, config)


def _load_target_manifest(
    manifest_path: Path,
    filters: Iterable[str],
) -> pd.DataFrame:
    manifest = read_csv_cached(str(manifest_path))
    if manifest is None:
        raise FileNotFoundError(f"Target manifest missing: {manifest_path}")
    if filters:
        manifest = manifest[manifest["Star Name"].isin(filters)]
    required_columns = {"Star Name", "Star Simbad Name"}
    missing = required_columns.difference(manifest.columns)
    if missing:
        raise ValueError(f"Target manifest missing required columns: {sorted(missing)}")
    return manifest.reset_index(drop=True)


def _build_base_payload(ephemeris, cadence: MinuteCadence) -> dict[str, np.ndarray]:
    saa_crossing = compute_saa_crossings(
        ephemeris.spacecraft_lat_deg, ephemeris.spacecraft_lon_deg
    )

    # Pre-compute datetime conversion once for all stars.
    # Store as datetime64[ns] (timestamp) so parquet writing is faster than
    # writing variable-length strings, and downstream code avoids re-parsing.
    mjd_array = np.asarray(cadence.mjd_utc, dtype=float)
    time_utc = Time(mjd_array, format="mjd", scale="utc")
    datetime_utc = time_utc.datetime64

    # SkyCoord objects for angular separation (kept for backward compat)
    earth_pc_sc = SkyCoord(
        ephemeris.earth_pc,
        unit=u.km,
        representation_type="cartesian",
    )
    sun_pc_sc = SkyCoord(
        ephemeris.sun_pc,
        unit=u.km,
        representation_type="cartesian",
    )
    moon_pc_sc = SkyCoord(
        ephemeris.moon_pc,
        unit=u.km,
        representation_type="cartesian",
    )

    # --- Precompute unit vectors for constraint evaluation (once for all stars) ---
    earth_pc_xyz = ephemeris.earth_pc  # (N, 3) km, s/c → Earth (nadir direction)
    sun_pc_xyz = ephemeris.sun_pc      # (N, 3) km, s/c → Sun
    moon_pc_xyz = ephemeris.moon_pc    # (N, 3) km, s/c → Moon

    nadir_unit = _normalise(earth_pc_xyz)           # s/c → Earth centre
    zenith_unit = -nadir_unit                       # Earth centre → s/c
    sun_unit = _normalise(sun_pc_xyz)
    moon_unit = _normalise(moon_pc_xyz)

    observer_dist_km = np.linalg.norm(earth_pc_xyz, axis=1)  # (N,)
    limb_angle_rad = np.arccos(
        np.clip(_R_EARTH_KM / observer_dist_km, -1.0, 1.0)
    )  # (N,)

    # Orbit boundary detection from sub-satellite latitude
    orbit_boundaries = detect_orbit_boundaries(ephemeris.spacecraft_lat_deg)
    orbit_slices = orbit_slices_from_boundaries(orbit_boundaries, len(mjd_array))

    return {
        "Time(MJD_UTC)": mjd_array,
        "Time_UTC": datetime_utc,
        "SAA_Crossing": np.round(saa_crossing, 1),
        # SkyCoord objects for backward-compatible separations
        "earth_pc": earth_pc_sc,
        "sun_pc": sun_pc_sc,
        "moon_pc": moon_pc_sc,
        # Unit vectors for constraint engine
        "nadir_unit": nadir_unit,
        "zenith_unit": zenith_unit,
        "sun_unit": sun_unit,
        "moon_unit": moon_unit,
        "observer_dist_km": observer_dist_km,
        "limb_angle_rad": limb_angle_rad,
        "orbit_slices": orbit_slices,
    }


def _build_star_visibility(
    payload: dict[str, np.ndarray],
    star_coord: SkyCoord,
    config: PandoraSchedulerConfig,
) -> pd.DataFrame:
    # --- Target unit vector (direction from observer to star) ---
    # Use the SkyCoord Earth-centre separation as the baseline Earth check.
    earth_center_sep_deg = payload["earth_pc"].separation(star_coord).deg

    # Build target unit vector in the same ECI frame used by the ephemeris.
    # SkyCoord.cartesian gives unit direction; replicate for each timestep.
    tgt_cart = star_coord.icrs.cartesian
    tgt_unit_1 = np.array([tgt_cart.x.value, tgt_cart.y.value, tgt_cart.z.value])
    tgt_unit_1 = tgt_unit_1 / np.linalg.norm(tgt_unit_1)
    N = len(payload["Time(MJD_UTC)"])
    target_unit = np.broadcast_to(tgt_unit_1, (N, 3)).copy()

    results = compute_visibility_with_constraints(
        target_unit=target_unit,
        nadir_unit=payload["nadir_unit"],
        sun_unit=payload["sun_unit"],
        moon_unit=payload["moon_unit"],
        observer_dist_km=payload["observer_dist_km"],
        zenith_unit=payload["zenith_unit"],
        limb_angle_rad=payload["limb_angle_rad"],
        orbit_slices=payload["orbit_slices"],
        earth_center_sep_deg=earth_center_sep_deg,
        config=config,
    )

    visible = results["visible"].astype(float)

    # Use pre-computed datetime array from payload (computed once for all stars)
    data = {
        "Time(MJD_UTC)": payload["Time(MJD_UTC)"],
        "Time_UTC": payload["Time_UTC"],
        "SAA_Crossing": payload["SAA_Crossing"],
        "Visible": np.round(visible, 1),
        "Earth_Sep": np.round(earth_center_sep_deg, 3),
        "Moon_Sep": np.round(results["moon_sep"], 3),
        "Sun_Sep": np.round(results["sun_sep"], 3),
        "Roll_Deg": np.round(results["roll_deg"], 2),
        "N_ST_Pass": results["n_st_pass"].astype(int),
    }
    return pd.DataFrame(data)


def _resolve_star_coord(
    row: pd.Series,
    star_metadata: dict[str, tuple[float, float]],
) -> SkyCoord:
    """Resolve star coordinates from catalog data only (no Simbad lookups)."""
    star_name = str(row.get("Star Name", ""))

    ra_val = row.get("RA")
    dec_val = row.get("DEC")

    # Use star_metadata as fallback if RA/DEC missing
    if (pd.isna(ra_val) or pd.isna(dec_val)) and star_name in star_metadata:
        fallback_ra, fallback_dec = star_metadata[star_name]
        if pd.isna(ra_val):
            ra_val = fallback_ra
        if pd.isna(dec_val):
            dec_val = fallback_dec

    if pd.notna(ra_val) and pd.notna(dec_val):
        return SkyCoord(
            ra=float(ra_val) * u.deg, dec=float(dec_val) * u.deg, frame="icrs"
        )

    # No Simbad fallback - raise error if coordinates not in catalog
    raise RuntimeError(f"No coordinates found in catalog for {star_name}")


def _build_star_metadata(manifest: pd.DataFrame) -> dict[str, tuple[float, float]]:
    if "RA" not in manifest.columns or "DEC" not in manifest.columns:
        return {}

    metadata: dict[str, tuple[float, float]] = {}
    for _, row in manifest.iterrows():
        ra = row.get("RA")
        dec = row.get("DEC")
        if pd.notna(ra) and pd.notna(dec):
            star_name = str(row.get("Star Name", ""))
            metadata[star_name] = (float(ra), float(dec))
    return metadata


def _build_planet_transits(
    manifest: pd.DataFrame,
    manifest_path: Path,
    output_root: Path,
    star_metadata: dict[str, tuple[float, float]],
    config: PandoraSchedulerConfig,
) -> list[tuple[str, str]]:
    if manifest.empty:
        return []

    required_columns = {
        "Planet Name",
        "Star Name",
        "Transit Duration (hrs)",
        "Period (days)",
        "Transit Epoch (BJD_TDB-2400000.5)",
    }
    missing = required_columns.difference(manifest.columns)
    if missing:
        LOGGER.info(
            "Manifest %s missing planet columns; skipping transit generation",
            manifest_path.name,
        )
        return []

    observer_location = EarthLocation(
        lat=0.0 * u.deg, lon=0.0 * u.deg, height=600.0 * u.km
    )

    generated: list[tuple[str, str]] = []
    skipped_existing: list[tuple[str, str]] = []

    for _, row in manifest.iterrows():
        star_name = str(row.get("Star Name", ""))
        planet_name = str(row.get("Planet Name", ""))

        star_visibility_path = (
            output_root / star_name / f"Visibility for {star_name}.parquet"
        )
        if not star_visibility_path.exists():
            LOGGER.warning(
                "Star visibility missing for %s; skipping planet %s",
                star_name,
                planet_name,
            )
            continue

        planet_dir = output_root / star_name / planet_name
        planet_dir.mkdir(parents=True, exist_ok=True)
        planet_output = planet_dir / f"Visibility for {planet_name}.parquet"
        if planet_output.exists() and not config.force_regenerate:
            skipped_existing.append((star_name, planet_name))
            generated.append((star_name, planet_name))
            continue

        planet_df = _compute_planet_transits(
            star_visibility_path,
            row,
            star_metadata,
            observer_location,
        )
        _write_visibility_parquet(planet_df, planet_output, config)
        if not planet_df.empty:
            generated.append((star_name, planet_name))

    if skipped_existing:
        example_pairs = ", ".join(
            f"{star}/{planet}" for star, planet in skipped_existing[:5]
        )
        LOGGER.info(
            "Skipping %d planet visibility files that already exist%s",
            len(skipped_existing),
            f" (e.g. {example_pairs})" if example_pairs else "",
        )

    return generated


def _compute_planet_transits(
    star_visibility_path: Path,
    planet_row: pd.Series,
    star_metadata: dict[str, tuple[float, float]],
    observer_location: EarthLocation,
) -> pd.DataFrame:
    star_visibility = read_parquet_cached(
        str(star_visibility_path),
        columns=["Time(MJD_UTC)", "Visible", "SAA_Crossing"],
    )
    if star_visibility is None or star_visibility.empty:
        raise FileNotFoundError(
            f"Star visibility missing or empty for {star_visibility_path}"
        )
    t_mjd = star_visibility["Time(MJD_UTC)"].to_numpy(dtype=float)
    visible_mask = star_visibility["Visible"].to_numpy(dtype=float)

    if t_mjd.size == 0:
        return pd.DataFrame(
            {
                "Transits": np.array([], dtype=float),
                "Transit_Start": np.array([], dtype=float),
                "Transit_Stop": np.array([], dtype=float),
                "Transit_Start_UTC": pd.Series([], dtype="datetime64[ns]"),
                "Transit_Stop_UTC": pd.Series([], dtype="datetime64[ns]"),
                "Transit_Coverage": np.array([], dtype=float),
                "SAA_Overlap": np.array([], dtype=float),
            }
        )

    transit_duration = planet_row["Transit Duration (hrs)"]
    period_days = planet_row["Period (days)"]
    epoch_bjd_tdb = planet_row["Transit Epoch (BJD_TDB-2400000.5)"]
    planet_name = planet_row["Planet Name"]

    if np.isnan(transit_duration) or np.isnan(period_days) or np.isnan(epoch_bjd_tdb):
        LOGGER.warning(
            "Incomplete ephemeris for %s; skipping planet visibility",
            planet_name,
        )
        return pd.DataFrame(
            {
                "Transits": np.array([], dtype=float),
                "Transit_Start": np.array([], dtype=float),
                "Transit_Stop": np.array([], dtype=float),
                "Transit_Start_UTC": pd.Series([], dtype="datetime64[ns]"),
                "Transit_Stop_UTC": pd.Series([], dtype="datetime64[ns]"),
                "Transit_Coverage": np.array([], dtype=float),
                "SAA_Overlap": np.array([], dtype=float),
            }
        )

    transit_duration = float(transit_duration) * u.hour
    period = float(period_days) * u.day

    star_coord = _resolve_star_coord(
        planet_row,
        star_metadata,
    )

    bjd_tdb = Time(
        float(epoch_bjd_tdb) + 2400000.5,
        format="jd",
        scale="tdb",
        location=observer_location,
    )
    light_time = bjd_tdb.light_travel_time(
        star_coord, kind="barycentric", location=observer_location
    )
    jd_tdb = bjd_tdb - light_time
    epoch_mjd_utc = Time(jd_tdb.mjd, format="mjd", scale="utc")

    half_obs_width = 0.75 * u.hour + np.maximum(
        1.0 * u.hour + transit_duration / 2.0, transit_duration
    )
    time_grid = Time(t_mjd, format="mjd", scale="utc")

    if period <= 0 * u.day:
        LOGGER.warning(
            "Non-positive period for %s; skipping planet visibility",
            planet_name,
        )
        return pd.DataFrame(
            {
                "Transits": np.array([], dtype=float),
                "Transit_Start": np.array([], dtype=float),
                "Transit_Stop": np.array([], dtype=float),
                "Transit_Start_UTC": pd.Series([], dtype="datetime64[ns]"),
                "Transit_Stop_UTC": pd.Series([], dtype="datetime64[ns]"),
                "Transit_Coverage": np.array([], dtype=float),
                "SAA_Overlap": np.array([], dtype=float),
            }
        )

    min_start_epoch = epoch_mjd_utc - half_obs_width
    elapsed_days = (time_grid[0] - min_start_epoch).to(u.day)
    min_pers_start = np.ceil((elapsed_days / period).value)

    first_transit = epoch_mjd_utc + min_pers_start * period

    mid_transits_list: list[Time] = []
    current = first_transit
    last_time = time_grid[-1]
    while current < last_time:
        mid_transits_list.append(current)
        current = current + period

    if not mid_transits_list:
        return pd.DataFrame(
            {
                "Transits": np.array([], dtype=float),
                "Transit_Start": np.array([], dtype=float),
                "Transit_Stop": np.array([], dtype=float),
                "Transit_Start_UTC": pd.Series([], dtype="datetime64[ns]"),
                "Transit_Stop_UTC": pd.Series([], dtype="datetime64[ns]"),
                "Transit_Coverage": np.array([], dtype=float),
                "SAA_Overlap": np.array([], dtype=float),
            }
        )

    mid_transits = Time(mid_transits_list)
    start_transits = mid_transits - transit_duration / 2.0
    end_transits = mid_transits + transit_duration / 2.0

    start_datetimes = start_transits.to_value("datetime")
    end_datetimes = end_transits.to_value("datetime")

    # Floor to minute precision using pandas vectorized operations
    start_datetimes = pd.to_datetime(start_datetimes).floor("min").to_pydatetime()
    end_datetimes = pd.to_datetime(end_datetimes).floor("min").to_pydatetime()

    saa_mask = star_visibility["SAA_Crossing"].to_numpy(dtype=float)
    T_mjd_utc = Time(t_mjd, format="mjd", scale="utc")
    T_iso_utc = Time(T_mjd_utc.iso, format="iso", scale="utc")
    dt_iso_utc = T_iso_utc.to_value("datetime")

    # Use boolean indexing for better performance
    dt_vis_times = dt_iso_utc[visible_mask == 1.0]
    dt_saa_times = dt_iso_utc[saa_mask == 1.0]

    coverage = np.zeros(len(start_datetimes), dtype=float)
    saa_overlap = np.zeros(len(start_datetimes), dtype=float)

    for idx, (start_dt, end_dt) in enumerate(zip(start_datetimes, end_datetimes)):
        # Use periods= to generate a half-open [start, end) range so the
        # denominator equals the number of 1-minute intervals, not N+1 (the
        # inclusive pd.date_range default would give off-by-one coverage).
        n_minutes = int((end_dt - start_dt).total_seconds() // 60)
        tran_minutes = pd.date_range(start_dt, periods=n_minutes, freq="min").to_pydatetime()
        if len(tran_minutes) == 0:
            continue
        minute_set = set(tran_minutes)
        tran_vis = minute_set.intersection(dt_vis_times)
        if len(tran_vis) > 0:
            coverage[idx] = len(tran_vis) / len(tran_minutes)

        saa_candidates = []
        for dt_val in dt_saa_times:
            if tran_minutes[0] <= dt_val <= tran_minutes[-1]:
                saa_candidates.append(dt_val)
        if saa_candidates:
            overlap = set(saa_candidates).intersection(tran_minutes)
            if overlap:
                saa_overlap[idx] = len(overlap) / len(tran_minutes)

    transit_df = pd.DataFrame(
        {
            "Transits": np.arange(len(start_datetimes), dtype=int),
            "Transit_Start": start_transits.value,
            "Transit_Stop": end_transits.value,
            "Transit_Start_UTC": start_datetimes,
            "Transit_Stop_UTC": end_datetimes,
            "Transit_Coverage": coverage,
            "SAA_Overlap": saa_overlap,
        }
    )
    return transit_df


def _apply_transit_overlaps(
    generated_planets: Iterable[tuple[str, str]],
    output_root: Path,
    config: PandoraSchedulerConfig | None = None,
) -> None:
    star_planets: dict[str, list[str]] = {}
    for star_name, planet_name in generated_planets:
        star_planets.setdefault(star_name, []).append(planet_name)

    for star_name, planets in star_planets.items():
        if len(planets) < 2:
            continue

        # Quick check: if all planet files already have Transit_Overlap, skip expensive recomputation
        all_have_overlap = True
        for planet in planets:
            planet_path = (
                output_root / star_name / planet / f"Visibility for {planet}.parquet"
            )
            if not planet_path.exists():
                all_have_overlap = False
                break
            # Quick check via parquet schema metadata.
            try:
                import pyarrow.parquet as pq

                schema = pq.read_schema(planet_path)
                if "Transit_Overlap" not in schema.names:
                    all_have_overlap = False
                    break
            except Exception:
                all_have_overlap = False
                break

        if all_have_overlap:
            continue  # Skip this star system - all planets already have overlaps computed

        planet_data: dict[str, pd.DataFrame] = {}
        minute_sets: dict[str, list[tuple[set, int]]] = {}

        for planet in planets:
            planet_path = (
                output_root / star_name / planet / f"Visibility for {planet}.parquet"
            )
            if not planet_path.exists():
                raise FileNotFoundError(
                    f"Expected planet visibility missing: {planet_path}"
                )
            df = read_parquet_cached(
                str(planet_path),
            )
            if df is None:
                raise FileNotFoundError(
                    f"Unable to read planet visibility: {planet_path}"
                )
            # Skip planets with no transits (empty DataFrame except for headers)
            if df.empty:
                continue
            planet_data[planet] = df
            sets: list[tuple[set, int]] = []

            # Vectorized datetime processing (much faster than iterrows)
            if (
                "Transit_Start_UTC" in df.columns
                and df["Transit_Start_UTC"].notna().any()
            ):
                # Use pre-existing datetime columns (fastest path)
                start_times = pd.to_datetime(df["Transit_Start_UTC"]).dt.floor("min")
                end_times = pd.to_datetime(df["Transit_Stop_UTC"]).dt.floor("min")
            else:
                # Fallback to MJD conversion (vectorized)
                start_mjd = df["Transit_Start"].to_numpy(dtype=float)
                end_mjd = df["Transit_Stop"].to_numpy(dtype=float)
                start_times = pd.Series(
                    Time(start_mjd, format="mjd", scale="utc").to_datetime()
                ).dt.floor("min")
                end_times = pd.Series(
                    Time(end_mjd, format="mjd", scale="utc").to_datetime()
                ).dt.floor("min")

            # Build minute sets for each transit (half-open [start, end)).
            for start_dt, end_dt in zip(start_times, end_times):
                n_minutes = int((end_dt - start_dt).total_seconds() // 60)
                minutes = list(
                    pd.date_range(start_dt, periods=n_minutes, freq="min").to_pydatetime()
                )
                if not minutes:
                    sets.append((set(), 0))
                else:
                    sets.append((set(minutes), len(minutes)))
            minute_sets[planet] = sets

        for planet, df in planet_data.items():
            overlaps = np.zeros(len(df), dtype=float)
            current_sets = minute_sets[planet]
            for idx, (minutes, total) in enumerate(current_sets):
                if total == 0:
                    continue
                best_overlap = 0.0
                for other_planet, other_sets in minute_sets.items():
                    if other_planet == planet:
                        continue
                    for other_minutes, other_total in other_sets:
                        if other_total == 0 or not other_minutes:
                            continue
                        shared = minutes.intersection(other_minutes)
                        if shared:
                            overlap_fraction = len(shared) / total
                            best_overlap = max(best_overlap, min(overlap_fraction, 1.0))
                overlaps[idx] = min(
                    best_overlap, 1.0
                )  # Ensure overlap never exceeds 1.0

            if "Transit_Overlap" in df.columns:
                df["Transit_Overlap"] = overlaps
            else:
                df["Transit_Overlap"] = overlaps

            planet_path = (
                output_root / star_name / planet / f"Visibility for {planet}.parquet"
            )
            if config is not None:
                _write_visibility_parquet(df, planet_path, config)
            else:
                df.to_parquet(
                    planet_path,
                    index=False,
                    engine="pyarrow",
                    compression="snappy",
                    write_statistics=False,
                    use_dictionary=False,
                )
