from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from Binary_Pulsar_Distances.atnf import read_atnf_long_with_errors
from Binary_Pulsar_Distances.eliminating import (
    add_dm_distance,
    add_gaia_distance,
    confirm_proper_motion,
    filter_binary,
    filter_in_globular,
    filter_position_uncertainty,
    matching_pipeline,
)

FIXTURES = Path(__file__).resolve().parent / "Binary_Pulsar_Distances" / "text_files_test"


class TestReadAtnfLongWithErrors:

    def test_parses_positions_pm_and_binary_type(self):
        """
        e2_input.csv has 4 pulsars in ATNF's "long with errors" format; check
        that positions, proper motions, and binary type are parsed correctly
        for a fully-populated row.
        """
        df = read_atnf_long_with_errors(FIXTURES / "e2_input.csv")

        assert list(df["jname"]) == ["J0002+6216", "J0006+1834", "J0011+08", "J0023+0923"]

        row = df[df["jname"] == "J0023+0923"].iloc[0]
        assert row["pmra_masyr"] == pytest.approx(-12.44)
        assert row["pmdec_masyr"] == pytest.approx(-6.16)
        assert row["binary_type"] == "ELL1"

        isolated = df[df["jname"] == "J0002+6216"].iloc[0]
        assert np.isnan(isolated["pmra_masyr"])
        assert isolated["binary_type"] == "*"

    def test_raj_err_converted_to_arcsec(self):
        """
        ATNF's RAJ_ERR is given in seconds of time; it should come out
        converted to arcsec (via 15 * cos(dec)), not left in time units.
        """
        df = read_atnf_long_with_errors(FIXTURES / "e2_input.csv")
        row = df[df["jname"] == "J0023+0923"].iloc[0]

        raj_err_time_s = 7.0e-06
        expected_arcsec = raj_err_time_s * 15.0 * np.cos(np.radians(row["decj_deg"]))
        assert row["raj_err_arcsec"] == pytest.approx(expected_arcsec)


class TestFilterPositionUncertainty:

    def test_keeps_only_sub_arcsec_on_both_axes(self):
        """
        Of the 4 pulsars in e2_input.csv, only 2 have both RA and Dec
        uncertainty under 1 arcsec.
        """
        df = read_atnf_long_with_errors(FIXTURES / "e2_input.csv")
        filtered = filter_position_uncertainty(df, max_arcsec=1.0)

        assert set(filtered["jname"]) == {"J0002+6216", "J0023+0923"}

    def test_excludes_unmeasured_uncertainty(self):
        """
        A pulsar with no measured position uncertainty (NaN) should be
        excluded, not treated as passing the cut.
        """
        df = pd.DataFrame(
            {
                "jname": ["J0000+0000"],
                "raj_err_arcsec": [np.nan],
                "decj_err_arcsec": [0.1],
            }
        )
        filtered = filter_position_uncertainty(df, max_arcsec=1.0)
        assert len(filtered) == 0


class TestFilterBinary:

    def test_keeps_only_binaries(self):
        df = read_atnf_long_with_errors(FIXTURES / "e2_input.csv")
        filtered = filter_binary(df)
        assert list(filtered["jname"]) == ["J0023+0923"]


class TestFilterInGlobular:

    def test_removes_known_globular_cluster_pulsar(self):
        """
        globular_test_long.csv has one pulsar in a globular cluster
        (J2140-2310B, a real entry in gc_pulsar_names.csv) and one that
        isn't (J0023+0923); only the latter should remain.
        """
        df = read_atnf_long_with_errors(FIXTURES / "globular_test_long.csv")
        filtered = filter_in_globular(df)
        assert list(filtered["jname"]) == ["J0023+0923"]


class TestConfirmProperMotion:

    def _synthetic_matches(self):
        return pd.DataFrame(
            {
                "jname": ["J0000+0000", "J0000+0000"],
                "source_id": [111, 222],
                "pmra_masyr": [10.0, 10.0],
                "pmra_err_masyr": [0.1, 0.1],
                "pmdec_masyr": [-5.0, -5.0],
                "pmdec_err_masyr": [0.1, 0.1],
                "pmra": [10.05, 20.0],
                "pmra_error": [0.1, 0.1],
                "pmdec": [-5.05, 5.0],
                "pmdec_error": [0.1, 0.1],
            }
        )

    def test_flags_agreeing_and_disagreeing_proper_motion(self):
        result = confirm_proper_motion(self._synthetic_matches(), n_sigma=3.0)

        agreeing = result[result["source_id"] == 111].iloc[0]
        disagreeing = result[result["source_id"] == 222].iloc[0]

        assert agreeing["pm_match"]
        assert not disagreeing["pm_match"]


class TestAddGaiaDistance:

    def test_prefers_gspphot_and_falls_back_to_parallax_inverse(self):
        matches = pd.DataFrame(
            {
                "parallax": [2.0, 1.0],
                "parallax_error": [0.1, 0.5],
                "distance_gspphot": [np.nan, 900.0],
                "distance_gspphot_lower": [np.nan, 850.0],
                "distance_gspphot_upper": [np.nan, 950.0],
            }
        )
        result = add_gaia_distance(matches)

        no_gspphot, has_gspphot = result.iloc[0], result.iloc[1]
        assert no_gspphot["distance_method"] == "parallax_inverse"
        assert no_gspphot["gaia_distance_pc"] == pytest.approx(500.0)
        assert has_gspphot["distance_method"] == "gspphot"
        assert has_gspphot["gaia_distance_pc"] == pytest.approx(900.0)


class TestAddDmDistance:

    def test_matches_published_vlbi_distance_for_j2222_0137(self):
        """
        PSR J2222-0137's distance has been independently measured via VLBI
        parallax at 267.3 +/- 1.2 pc (Guo et al. 2021). The YMW16 DM-distance
        computed here should land in the same ballpark (DM-distance methods
        are typically only good to ~20%, so this isn't a tight tolerance).
        """
        pytest.importorskip("pygedm")

        matches = pd.DataFrame(
            {"gl_deg": [62.018572], "gb_deg": [-46.075389], "dm": [3.2826]}
        )
        result = add_dm_distance(matches)

        assert result["dm_distance_pc"].iloc[0] == pytest.approx(267.3, rel=0.2)


class TestMatchingPipeline:

    def test_pipeline_runs_end_to_end_against_a_real_pulsar(self, tmp_path):
        """
        Runs the full pipeline against a single real, well-known binary
        pulsar (J2222-0137, ATNF parameters fetched live via psrqpy) with
        real, live Gaia queries -- a plumbing/smoke test, not a check that
        this specific companion is Gaia-detectable (many pulsar white-dwarf
        companions are too faint for Gaia's sensitivity limit, so pm_match
        may legitimately be False for every candidate here).
        """
        try:
            import pygedm  # noqa: F401

            include_dm_distance = True
        except ImportError:
            include_dm_distance = False

        input_file = FIXTURES / "small_atnf_known_binary.csv"
        output_file = tmp_path / "matches.csv"

        matches = matching_pipeline(
            str(input_file),
            str(output_file),
            radius_arcsec=60.0,
            include_dm_distance=include_dm_distance,
        )

        assert output_file.exists()
        if len(matches) > 0:
            expected_columns = ["pm_match", "pm_sigma_ra", "gaia_distance_pc", "distance_method"]
            if include_dm_distance:
                expected_columns.append("dm_distance_pc")
            for column in expected_columns:
                assert column in matches.columns
