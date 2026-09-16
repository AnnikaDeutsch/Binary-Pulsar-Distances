"""Parsing for ATNF psrcat catalogue exports."""

import numpy as np
import pandas as pd
from astropy.coordinates import Angle, SkyCoord
import astropy.units as u

# Column positions (after whitespace-splitting one line) in ATNF psrcat's
# "long with errors" printout: each parameter is followed by its
# uncertainty and a reference code, e.g. "RAJ RAJ_ERR RAJ_REF".
_FIELD_INDEX = {
    "jname": 1,
    "raj": 3,
    "raj_err": 4,
    "decj": 6,
    "decj_err": 7,
    "pmra": 9,
    "pmra_err": 10,
    "pmdec": 12,
    "pmdec_err": 13,
    "posepoch": 15,
    "dm": 17,
    "dm_err": 18,
    "binary_type": 20,
}


def _parse_value(raw):
    """Parse an ATNF numeric field, treating '*' (unmeasured) and '0'
    (ATNF's convention for an unmeasured uncertainty) as missing."""
    return np.nan if raw in ("*", "0") else float(raw)


def read_atnf_long_with_errors(path):
    """Parse an ATNF psrcat "long with errors" export into a DataFrame.

    Each pulsar becomes one row with columns: jname, raj_deg, raj_err_arcsec,
    decj_deg, decj_err_arcsec, pmra_masyr, pmra_err_masyr, pmdec_masyr,
    pmdec_err_masyr, posepoch_mjd, dm, dm_err, binary_type, gl_deg, gb_deg.

    RAJ_ERR is given by ATNF in seconds of time (since RAJ's sexagesimal
    unit is hours:minutes:seconds of time); it's converted here to arcsec
    (via 15 * cos(dec)) so it's directly comparable to DECJ_ERR, which ATNF
    already gives in arcsec.

    Args:
        path (str): Path to the ATNF "long with errors" text export.

    Returns:
        pandas.DataFrame: One row per pulsar.
    """
    rows = []
    with open(path) as f:
        for line in f:
            fields = line.split()
            if not fields:
                continue

            raj_angle = Angle(fields[_FIELD_INDEX["raj"]] + " hours")
            decj_angle = Angle(fields[_FIELD_INDEX["decj"]] + " degrees")

            raj_err_time_s = _parse_value(fields[_FIELD_INDEX["raj_err"]])
            raj_err_arcsec = (
                raj_err_time_s * 15.0 * np.cos(decj_angle.radian)
                if not np.isnan(raj_err_time_s)
                else np.nan
            )

            rows.append(
                {
                    "jname": fields[_FIELD_INDEX["jname"]],
                    "raj_deg": raj_angle.degree,
                    "raj_err_arcsec": raj_err_arcsec,
                    "decj_deg": decj_angle.degree,
                    "decj_err_arcsec": _parse_value(fields[_FIELD_INDEX["decj_err"]]),
                    "pmra_masyr": _parse_value(fields[_FIELD_INDEX["pmra"]]),
                    "pmra_err_masyr": _parse_value(fields[_FIELD_INDEX["pmra_err"]]),
                    "pmdec_masyr": _parse_value(fields[_FIELD_INDEX["pmdec"]]),
                    "pmdec_err_masyr": _parse_value(fields[_FIELD_INDEX["pmdec_err"]]),
                    "posepoch_mjd": _parse_value(fields[_FIELD_INDEX["posepoch"]]),
                    "dm": _parse_value(fields[_FIELD_INDEX["dm"]]),
                    "dm_err": _parse_value(fields[_FIELD_INDEX["dm_err"]]),
                    "binary_type": fields[_FIELD_INDEX["binary_type"]],
                }
            )

    df = pd.DataFrame(rows)
    galactic = SkyCoord(
        ra=df["raj_deg"].to_numpy() * u.deg,
        dec=df["decj_deg"].to_numpy() * u.deg,
        frame="icrs",
    ).galactic
    df["gl_deg"] = galactic.l.degree
    df["gb_deg"] = galactic.b.degree
    return df
