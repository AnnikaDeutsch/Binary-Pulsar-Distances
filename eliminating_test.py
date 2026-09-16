from pathlib import Path

import pytest

from Binary_Pulsar_Distances.eliminating import (
    check_binary,
    check_in_globular,
    check_pos_uncertainty,
    matching_pipeline,
)

FIXTURES = Path(__file__).resolve().parent / "Binary_Pulsar_Distances" / "text_files_test"


def _count_lines(path):
    with open(path) as f:
        return sum(1 for _ in f)


class TestBinaryCheck:

    @pytest.mark.xfail(
        reason=(
            "check_binary() indexes into the check_pos_uncertainty()-output "
            "column layout (binary flag at index 11), but this fixture uses "
            "the plain ATNF export layout (binary flag at index 8), so it "
            "raises IndexError. Needs a decision on which layout check_binary "
            "should support (see project discussion) before this is fixed."
        ),
        strict=True,
    )
    def test_eliminates_all_not_binaries(self, tmp_path):
        """
        Given a small file with a known number of isolated pulsars, check_binary()
        should eliminate all of them and keep only the binary pulsars.
        """
        input_file = FIXTURES / "e1_input.csv"
        output_file = tmp_path / "e1_output.csv"
        check_binary(str(input_file), str(output_file))

        assert _count_lines(output_file) == 2


class TestPosUncertaintyCheck:

    def test_pos_uncertainty_check(self, tmp_path):
        """
        Given a small file (e2_input.csv) with 2 acceptable and 2 unacceptable
        position uncertainties, the output file should have only 2 entries.
        """
        input_file = FIXTURES / "e2_input.csv"
        output_file = tmp_path / "e2_output.csv"
        check_pos_uncertainty(str(input_file), str(output_file))

        assert _count_lines(output_file) == 2


class TestInGlobularCheck:

    def test_in_globular_check(self, tmp_path):
        """
        Given a file with 2 pulsars not in globular clusters and 1 pulsar in a
        globular cluster, the globular-cluster pulsar should be removed, leaving
        a file with the 2 remaining pulsars.
        """
        input_file = FIXTURES / "e3_input.csv"
        output_file = tmp_path / "e3_output.csv"
        check_in_globular(str(input_file), str(output_file))

        assert _count_lines(output_file) == 2


class TestMatchingPipeline:

    @pytest.mark.skip(
        reason=(
            "Runs the full pipeline, including a live Gaia cone search, "
            "against the entire all_atnf.csv fixture (~4000 pulsars) -- far "
            "too slow for a routine test run (observed to take well over 2 "
            "minutes without finishing). Needs a small hand-picked fixture "
            "(a handful of pulsars) before this can run routinely; also "
            "blocked on the check_binary() column-layout bug (see "
            "TestBinaryCheck)."
        ),
    )
    def test_on_all_in_atnf_with_dr3(self, tmp_path):
        """
        The full matching_pipeline() should return the same number of Gaia
        matches as the reference output file.
        """
        input_file = FIXTURES / "all_atnf.csv"
        output_file = tmp_path / "t1_output.csv"
        no_pos = tmp_path / "t1_nopos.csv"
        no_bin = tmp_path / "t1_no_bin.csv"
        no_glob = tmp_path / "t1_noglob.csv"
        match_all_params = tmp_path / "t1_matchall.csv"

        matching_pipeline(
            str(input_file),
            str(output_file),
            str(no_pos),
            str(no_bin),
            str(no_glob),
            str(match_all_params),
        )

        expected_file = FIXTURES / "all_final_short.csv"
        assert _count_lines(output_file) == _count_lines(expected_file)
