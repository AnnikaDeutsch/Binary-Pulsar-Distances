"""PSRmatch: cross-match pulsar catalogues against optical surveys."""

from .atnf import read_atnf_long_with_errors
from .eliminating import (
    add_dm_distance,
    add_gaia_distance,
    confirm_proper_motion,
    filter_binary,
    filter_has_proper_motion,
    filter_in_globular,
    filter_position_uncertainty,
    get_matches,
    matching_pipeline,
    pretty_print,
    pretty_print_matches,
    psr_to_gaia,
)

__all__ = [
    "add_dm_distance",
    "add_gaia_distance",
    "confirm_proper_motion",
    "filter_binary",
    "filter_has_proper_motion",
    "filter_in_globular",
    "filter_position_uncertainty",
    "get_matches",
    "matching_pipeline",
    "pretty_print",
    "pretty_print_matches",
    "psr_to_gaia",
    "read_atnf_long_with_errors",
]
