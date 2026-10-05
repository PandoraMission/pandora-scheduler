"""Tests for xml/parameters.py VDA branch coverage."""

import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd
import pytest

from pandorascheduler_rework.xml.parameters import (
    populate_nirda_parameters,
    populate_vda_parameters,
)


# ---------------------------------------------------------------------------
# populate_nirda_parameters
# ---------------------------------------------------------------------------


class TestPopulateNirdaEdgeCases:
    """NIRDA branches not covered by test_xml.py."""

    def test_empty_targ_info_writes_empty_element(self):
        """Empty targ_info emits <AcquireInfCamImages/> with a warning."""
        root = ET.Element("PayloadParameters")
        populate_nirda_parameters(root, pd.DataFrame(), 3600.0)
        child = root.find("AcquireInfCamImages")
        assert child is not None
        assert len(child) == 0

    def test_no_nirda_columns_skips(self):
        """DataFrame with no NIRDA_* columns skips population."""
        root = ET.Element("PayloadParameters")
        targ_info = pd.DataFrame({"Star Name": ["TestStar"], "RA": [0.0]})
        populate_nirda_parameters(root, targ_info, 3600.0)
        elem = root.find("AcquireInfCamImages")
        assert elem is not None
        assert len(elem) == 0


# ---------------------------------------------------------------------------
# populate_vda_parameters — branch coverage
# ---------------------------------------------------------------------------


class TestPopulateVdaTargetIdRaDec:
    """VDA TargetID, TargetRA, TargetDEC branches."""

    @staticmethod
    def _base_vda_row(**overrides):
        base = {
            "Star Name": "TestStar",
            "Planet Name": "TestStar b",
            "RA": 120.0,
            "DEC": 30.0,
            "VDA_TargetID": "SET_BY_SCHEDULER",
            "VDA_TargetRA": 0.0,
            "VDA_TargetDEC": 0.0,
            "StarRoiDetMethod": 2,
            "VDA_StarRoiDetMethod": 2,
        }
        base.update(overrides)
        return pd.DataFrame([base])

    def test_target_id_uses_identifier(self):
        root = ET.Element("PP")
        targ = self._base_vda_row()
        populate_vda_parameters(root, targ, 3600.0)
        vda = root.find("AcquireVisCamScienceData")
        tid = vda.find("TargetID")
        assert tid is not None
        assert tid.text  # non-empty identifier

    def test_target_ra_dec_populated(self):
        root = ET.Element("PP")
        targ = self._base_vda_row()
        populate_vda_parameters(root, targ, 3600.0)
        vda = root.find("AcquireVisCamScienceData")
        assert vda.find("TargetRA") is not None
        assert vda.find("TargetDEC") is not None

    def test_missing_ra_warns(self, caplog):
        root = ET.Element("PP")
        targ = self._base_vda_row(RA=float("nan"))
        populate_vda_parameters(root, targ, 3600.0)
        assert "Missing RA" in caplog.text

    def test_missing_dec_warns(self, caplog):
        root = ET.Element("PP")
        targ = self._base_vda_row(DEC=float("nan"))
        populate_vda_parameters(root, targ, 3600.0)
        assert "Missing DEC" in caplog.text


class TestPopulateVdaStarRoiDetMethod:
    """VDA StarRoiDetMethod branch coverage."""

    @staticmethod
    def _row(**overrides):
        base = {
            "Star Name": "X",
            "Planet Name": "X b",
            "RA": 0.0,
            "DEC": 0.0,
            "VDA_StarRoiDetMethod": 2,
            "StarRoiDetMethod": 2,
        }
        base.update(overrides)
        return pd.DataFrame([base])

    def test_method_2(self):
        root = ET.Element("PP")
        populate_vda_parameters(root, self._row(), 3600.0)
        vda = root.find("AcquireVisCamScienceData")
        method = vda.find("StarRoiDetMethod")
        assert method is not None
        assert method.text == "2"

    def test_method_placeholder_falls_back(self):
        root = ET.Element("PP")
        targ = self._row(VDA_StarRoiDetMethod="SET_BY_TARGET_DEFINITION_FILE", StarRoiDetMethod=1)
        populate_vda_parameters(root, targ, 3600.0)
        vda = root.find("AcquireVisCamScienceData")
        method = vda.find("StarRoiDetMethod")
        assert method is not None
        assert method.text == "1"

    def test_method_nan_raises(self):
        """NaN fallback StarRoiDetMethod raises when VDA column is a placeholder."""
        root = ET.Element("PP")
        targ = self._row(
            VDA_StarRoiDetMethod="SET_BY_TARGET_DEFINITION_FILE",
            StarRoiDetMethod=float("nan"),
        )
        with pytest.raises(ValueError, match="NaN StarRoiDetMethod"):
            populate_vda_parameters(root, targ, 3600.0)

    def test_method_unresolved_placeholder_raises(self):
        root = ET.Element("PP")
        targ = self._row(
            VDA_StarRoiDetMethod="SET_BY_TARGET_DEFINITION_FILE",
            StarRoiDetMethod="SET_BY_TARGET_DEFINITION_FILE",
        )
        with pytest.raises(ValueError, match="Unresolved placeholder"):
            populate_vda_parameters(root, targ, 3600.0)


class TestPopulateVdaMaxNumStarRois:
    """VDA MaxNumStarRois conditional logic."""

    @staticmethod
    def _row(**overrides):
        base = {
            "Star Name": "S",
            "Planet Name": "S b",
            "RA": 0.0,
            "DEC": 0.0,
            "VDA_MaxNumStarRois": 5,
            "VDA_StarRoiDetMethod": 1,
            "StarRoiDetMethod": 1,
            "numPredefinedStarRois": 3,
        }
        base.update(overrides)
        return pd.DataFrame([base])

    def test_method_1_uses_predefined(self):
        root = ET.Element("PP")
        populate_vda_parameters(root, self._row(), 3600.0)
        vda = root.find("AcquireVisCamScienceData")
        elem = vda.find("MaxNumStarRois")
        assert elem is not None
        assert elem.text == "3"  # numPredefinedStarRois value

    def test_method_2_uses_9(self):
        root = ET.Element("PP")
        populate_vda_parameters(root, self._row(StarRoiDetMethod=2, VDA_StarRoiDetMethod=2), 3600.0)
        vda = root.find("AcquireVisCamScienceData")
        elem = vda.find("MaxNumStarRois")
        assert elem is not None
        assert elem.text == "9"


class TestPopulateVdaNumPredefinedStarRois:
    """VDA numPredefinedStarRois branch."""

    @staticmethod
    def _row(**overrides):
        base = {
            "Star Name": "S",
            "Planet Name": "S b",
            "RA": 0.0,
            "DEC": 0.0,
            "VDA_numPredefinedStarRois": 0,
            "VDA_StarRoiDetMethod": 1,
            "StarRoiDetMethod": 1,
            "numPredefinedStarRois": 3,
        }
        base.update(overrides)
        return pd.DataFrame([base])

    def test_method_1_emits_value(self):
        root = ET.Element("PP")
        populate_vda_parameters(root, self._row(), 3600.0)
        vda = root.find("AcquireVisCamScienceData")
        elem = vda.find("numPredefinedStarRois")
        assert elem is not None
        assert elem.text == "3"

    def test_method_2_skips(self):
        root = ET.Element("PP")
        populate_vda_parameters(
            root, self._row(StarRoiDetMethod=2, VDA_StarRoiDetMethod=2), 3600.0
        )
        vda = root.find("AcquireVisCamScienceData")
        elem = vda.find("numPredefinedStarRois")
        assert elem is None


class TestPopulateVdaNumTotalFrames:
    """VDA NumTotalFramesRequested computation."""

    @staticmethod
    def _row(**overrides):
        base = {
            "Star Name": "S",
            "Planet Name": "S b",
            "RA": 0.0,
            "DEC": 0.0,
            "VDA_NumTotalFramesRequested": 0,
            "VDA_ExposureTime_us": 10000,  # 10ms = 0.01s
            "VDA_FramesPerCoadd": 10,
            "VDA_StarRoiDetMethod": 2,
            "StarRoiDetMethod": 2,
        }
        base.update(overrides)
        return pd.DataFrame([base])

    def test_frame_count_computed(self):
        root = ET.Element("PP")
        # 3600s / (0.01s * 10) = 36000 frames; floor / coadd * coadd = 36000
        populate_vda_parameters(root, self._row(), 3600.0)
        vda = root.find("AcquireVisCamScienceData")
        elem = vda.find("NumTotalFramesRequested")
        assert elem is not None
        frames = int(elem.text)
        assert frames > 0
        assert frames % 10 == 0  # multiple of FramesPerCoadd

    def test_zero_exposure_skips(self):
        root = ET.Element("PP")
        populate_vda_parameters(root, self._row(VDA_ExposureTime_us=0), 3600.0)
        vda = root.find("AcquireVisCamScienceData")
        elem = vda.find("NumTotalFramesRequested")
        # Zero exposure time → skip (doesn't write element)
        assert elem is None

    def test_nan_exposure_skips(self):
        root = ET.Element("PP")
        populate_vda_parameters(
            root, self._row(VDA_ExposureTime_us=float("nan")), 3600.0
        )
        vda = root.find("AcquireVisCamScienceData")
        elem = vda.find("NumTotalFramesRequested")
        assert elem is None


class TestPopulateVdaPredefinedRoiCoords:
    """VDA PredefinedStarRoiRa/Dec with ROI_coord_* columns."""

    @staticmethod
    def _row(**overrides):
        base = {
            "Star Name": "S",
            "Planet Name": "S b",
            "RA": 100.0,
            "DEC": -20.0,
            "VDA_PredefinedStarRoiRa": 0,
            "VDA_PredefinedStarRoiDec": 0,
            "VDA_StarRoiDetMethod": 1,
            "StarRoiDetMethod": 1,
            "numPredefinedStarRois": 2,
            "ROI_coord_1": "[100.1, -20.1]",
            "ROI_coord_2": "[101.0, -21.0]",
        }
        base.update(overrides)
        return pd.DataFrame([base])

    def test_roi_coords_emitted(self):
        root = ET.Element("PP")
        populate_vda_parameters(root, self._row(), 3600.0)
        vda = root.find("AcquireVisCamScienceData")
        ra_elem = vda.find("PredefinedStarRoiRa")
        dec_elem = vda.find("PredefinedStarRoiDec")
        assert ra_elem is not None
        assert dec_elem is not None
        # First slot should be overwritten with target RA/DEC
        ra1 = ra_elem.find("RA1")
        assert ra1 is not None
        assert float(ra1.text) == pytest.approx(100.0, abs=0.01)

    def test_method_2_skips_roi_coords(self):
        root = ET.Element("PP")
        populate_vda_parameters(
            root,
            self._row(StarRoiDetMethod=2, VDA_StarRoiDetMethod=2),
            3600.0,
        )
        vda = root.find("AcquireVisCamScienceData")
        assert vda.find("PredefinedStarRoiRa") is None
        assert vda.find("PredefinedStarRoiDec") is None

    def test_missing_roi_coord_columns_warns(self, tmp_path, caplog):
        root = ET.Element("PP")
        row_data = {
            "Star Name": "S",
            "Planet Name": "S b",
            "RA": 0.0,
            "DEC": 0.0,
            "VDA_PredefinedStarRoiRa": 0,
            "VDA_PredefinedStarRoiDec": 0,
            "VDA_StarRoiDetMethod": 1,
            "StarRoiDetMethod": 1,
        }
        populate_vda_parameters(root, pd.DataFrame([row_data]), 3600.0)
        assert "Missing ROI_coord" in caplog.text


class TestPopulateVdaEmpty:
    """VDA with empty targ_info."""

    def test_empty_writes_stub(self):
        root = ET.Element("PP")
        populate_vda_parameters(root, pd.DataFrame(), 3600.0)
        elem = root.find("AcquireVisCamScienceData")
        assert elem is not None
        assert len(elem) == 0
