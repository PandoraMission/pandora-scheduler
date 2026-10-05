"""Unit tests for the XML module."""

import xml.etree.ElementTree as ET
from datetime import datetime

import pandas as pd

from pandorascheduler_rework.xml.builder import (
    _build_observational_parameters,
    _duration_in_seconds,
    observation_sequence,
)
from pandorascheduler_rework.xml.parameters import (
    populate_nirda_parameters,
    populate_vda_parameters,
)


class TestObservationSequence:
    def test_creates_observation_sequence_element(self):
        visit = ET.Element("Visit")
        targ_info = pd.DataFrame()
        
        seq = observation_sequence(
            visit=visit,
            obs_seq_id="SEQ001",
            target_name="TestTarget",
            priority="High",
            start=datetime(2026, 1, 1, 12, 0, 0),
            stop=datetime(2026, 1, 1, 13, 0, 0),
            ra=10.5,
            dec=-20.5,
            targ_info=targ_info,
        )
        
        assert seq.tag == "Observation_Sequence"
        assert seq.find("ID").text == "SEQ001"
        assert seq.find("Observational_Parameters/Target").text == "TestTarget"

    def test_handles_string_timestamps(self):
        visit = ET.Element("Visit")
        targ_info = pd.DataFrame()
        
        seq = observation_sequence(
            visit=visit,
            obs_seq_id="SEQ002",
            target_name="TestTarget",
            priority="High",
            start="2026-01-01T12:00:00Z",
            stop="2026-01-01T13:00:00Z",
            ra=10.5,
            dec=-20.5,
            targ_info=targ_info,
        )
        
        timing = seq.find("Observational_Parameters/Timing")
        assert timing.find("Start").text == "2026-01-01T12:00:00Z"
        assert timing.find("Stop").text == "2026-01-01T13:00:00Z"


class TestHelperFunctions:
    def test_build_observational_parameters(self):
        params = _build_observational_parameters(
            "TargetA", "Low", 
            datetime(2026, 1, 1, 10, 0), 
            datetime(2026, 1, 1, 11, 0),
            100.0, 45.0
        )
        
        assert params["Target"] == "TargetA"
        assert params["Priority"] == "Low"
        assert params["Timing"][2] == "2026-01-01T10:00:00Z"
        assert params["Boresight"][2] == "100.0"

    def test_duration_in_seconds(self):
        start = datetime(2026, 1, 1, 10, 0, 0)
        stop = datetime(2026, 1, 1, 10, 1, 30)
        assert _duration_in_seconds(start, stop) == 90.0

    def test_duration_in_seconds_strings(self):
        start = "2026-01-01T10:00:00Z"
        stop = "2026-01-01T10:01:30Z"
        assert _duration_in_seconds(start, stop) == 90.0


class TestPopulateNirdaParameters:
    def test_populates_basic_nirda_parameters(self):
        root = ET.Element("Root")
        targ_info = pd.DataFrame([{
            "NIRDA_Filter": "Open",
            "NIRDA_ReadoutMode": "Fast",
            # NIRDA_IntegrationTime_s is in ignore list, so it won't be in XML
            "NIRDA_IntegrationTime_s": 0.5,
        }])
        
        populate_nirda_parameters(root, targ_info, diff_in_seconds=10.0)
        
        nirda = root.find("AcquireInfCamImages")
        assert nirda is not None
        assert nirda.find("Filter").text == "Open"
        assert nirda.find("ReadoutMode").text == "Fast"
        # IntegrationTime_s is intentionally ignored in XML generation
        assert nirda.find("IntegrationTime_s") is None

    def test_calculates_sc_integrations(self):
        root = ET.Element("Root")
        targ_info = pd.DataFrame([{
            "NIRDA_SC_Integrations": "CALCULATE",
            "NIRDA_IntegrationTime_s": 2.0,
        }])
        
        populate_nirda_parameters(root, targ_info, diff_in_seconds=10.0)
        
        nirda = root.find("AcquireInfCamImages")
        # 10s / 2s = 5 integrations
        assert nirda.find("SC_Integrations").text == "5"

    def test_handles_empty_dataframe(self):
        root = ET.Element("Root")
        targ_info = pd.DataFrame()
        
        populate_nirda_parameters(root, targ_info, diff_in_seconds=10.0)
        
        nirda = root.find("AcquireInfCamImages")
        assert nirda is not None
        assert len(list(nirda)) == 0  # Should be empty


class TestPopulateVdaParameters:
    def test_populates_basic_vda_parameters(self):
        root = ET.Element("Root")
        targ_info = pd.DataFrame([{
            "VDA_Filter": "Visible",
            "VDA_ExposureTime_us": 100000,
        }])
        
        populate_vda_parameters(root, targ_info, diff_in_seconds=10.0)
        
        vda = root.find("AcquireVisCamScienceData")
        assert vda is not None
        assert vda.find("Filter").text == "Visible"
        assert vda.find("ExposureTime_us").text == "100000"

    def test_uses_ra_dec_from_dataframe(self):
        root = ET.Element("Root")
        targ_info = pd.DataFrame([{
            "VDA_TargetRA": "USE_DATAFRAME",
            "VDA_TargetDEC": "USE_DATAFRAME",
            "RA": 123.45,
            "DEC": -67.89,
        }])
        
        populate_vda_parameters(root, targ_info, diff_in_seconds=10.0)
        
        vda = root.find("AcquireVisCamScienceData")
        assert vda.find("TargetRA").text == "123.45"
        assert vda.find("TargetDEC").text == "-67.89"

    def test_calculates_num_total_frames(self):
        root = ET.Element("Root")
        targ_info = pd.DataFrame([{
            "VDA_NumTotalFramesRequested": "CALCULATE",
            "VDA_ExposureTime_us": 500000,  # 0.5s
            "VDA_FramesPerCoadd": 10,
        }])
        
        # 10s duration / 0.5s exposure = 20 frames
        # 20 frames / 10 coadd = 2 coadds
        # 2 coadds * 10 frames = 20 frames total
        populate_vda_parameters(root, targ_info, diff_in_seconds=10.0)
        
        vda = root.find("AcquireVisCamScienceData")
        assert vda.find("NumTotalFramesRequested").text == "20"

    def test_handles_star_roi_det_method(self):
        root = ET.Element("Root")
        targ_info = pd.DataFrame([{
            "VDA_StarRoiDetMethod": "SET_BY_TARGET_DEFINITION_FILE",
            "StarRoiDetMethod": 2,
        }])
        
        populate_vda_parameters(root, targ_info, diff_in_seconds=10.0)
        
        vda = root.find("AcquireVisCamScienceData")
        assert vda.find("StarRoiDetMethod").text == "2"

    def test_predefined_star_roi_slot_zero_matches_target_ra_dec(self):
        root = ET.Element("Root")
        targ_info = pd.DataFrame([
            {
                "VDA_StarRoiDetMethod": "SET_BY_TARGET_DEFINITION_FILE",
                "StarRoiDetMethod": 1,
                "VDA_numPredefinedStarRois": "SET_BY_TARGET_DEFINITION_FILE",
                "numPredefinedStarRois": 2,
                "VDA_PredefinedStarRoiRa": "SET_BY_TARGET_DEFINITION_FILE",
                "VDA_PredefinedStarRoiDec": "SET_BY_TARGET_DEFINITION_FILE",
                "RA": 123.45,
                "DEC": -67.89,
                "ROI_coord_0": "[1.0, 2.0]",
                "ROI_coord_1": "[3.0, 4.0]",
            }
        ])

        populate_vda_parameters(root, targ_info, diff_in_seconds=10.0)

        vda = root.find("AcquireVisCamScienceData")
        assert vda is not None
        roi_ra = vda.find("PredefinedStarRoiRa")
        roi_dec = vda.find("PredefinedStarRoiDec")
        assert roi_ra is not None
        assert roi_dec is not None
        assert roi_ra.find("RA1").text == "123.450000"
        assert roi_dec.find("Dec1").text == "-67.890000"
        assert roi_ra.find("RA2").text == "3.000000"
        assert roi_dec.find("Dec2").text == "4.000000"


class TestXmlIntegration:
    def test_full_observation_sequence_with_all_parameters(self):
        visit = ET.Element("Visit")
        targ_info = pd.DataFrame([{
            "NIRDA_Filter": "Open",
            "NIRDA_TargetID": "USE_PLANET_NAME",
            "VDA_Filter": "Visible",
            "Planet Name": "WASP-121 b",
        }])
        
        seq = observation_sequence(
            visit=visit,
            obs_seq_id="SEQ_FULL",
            target_name="WASP-121b",
            priority="High",
            start=datetime(2026, 1, 1, 12, 0),
            stop=datetime(2026, 1, 1, 13, 0),
            ra=10.0,
            dec=-10.0,
            targ_info=targ_info,
        )
        
        # Check structure
        assert seq.find("ID").text == "SEQ_FULL"
        assert seq.find("Observational_Parameters/Target").text == "WASP-121b"
        
        # Check payload
        nirda = seq.find("Payload_Parameters/AcquireInfCamImages")
        assert nirda.find("Filter").text == "Open"
        assert nirda.find("TargetID").text == "WASP-121b"  # Space removed, suffix kept
        
        vda = seq.find("Payload_Parameters/AcquireVisCamScienceData")
        assert vda.find("Filter").text == "Visible"
