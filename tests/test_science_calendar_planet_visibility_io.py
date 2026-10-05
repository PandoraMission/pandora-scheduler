import pandas as pd
from astropy.time import Time

from pandorascheduler_rework.science_calendar import (
    _read_planet_visibility,
    _transit_windows,
)


def test_read_planet_visibility_parquet_uses_transit_columns(tmp_path, caplog):
    """Regression test: planet visibility parquet contains transit windows, not timeline.

    Historically, science_calendar tried to read planet parquet files with columns
    ["Time(MJD_UTC)", "Visible"], which fails for transit-window parquet schemas.
    """

    planet_name = "TOI-2076b"
    directory = tmp_path / "targets" / "TOI-2076" / planet_name
    directory.mkdir(parents=True)

    mjd_start = Time("2026-02-05T00:00:00.750", format="isot", scale="utc").mjd
    mjd_stop = Time("2026-02-05T01:00:00.250", format="isot", scale="utc").mjd

    df = pd.DataFrame({"Transit_Start": [mjd_start], "Transit_Stop": [mjd_stop]})
    parquet_path = directory / f"Visibility for {planet_name}.parquet"
    df.to_parquet(parquet_path)

    loaded = _read_planet_visibility(directory, planet_name)
    assert loaded is not None
    assert list(loaded.columns) == ["Transit_Start", "Transit_Stop"]

    windows = _transit_windows(loaded)
    assert windows is not None
    start_times, stop_times = windows

    # Verify _transit_windows still does the expected rounding
    assert start_times[0].isoformat(sep=" ") == "2026-02-05 00:00:01"
    assert stop_times[0].isoformat(sep=" ") == "2026-02-05 01:00:00"

    # Ensure we did not emit the pyarrow schema mismatch error.
    assert not any(
        "No match for FieldRef.Name(Time(MJD_UTC))" in record.getMessage()
        for record in caplog.records
    )
