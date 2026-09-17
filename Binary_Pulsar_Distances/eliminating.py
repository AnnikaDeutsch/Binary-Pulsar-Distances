import os
import time
from concurrent.futures import ThreadPoolExecutor, TimeoutError as FutureTimeoutError

import numpy as np
import pandas as pd
import requests
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astroquery.gaia import Gaia

from .atnf import read_atnf_long_with_errors

# astroquery/the Gaia TAP+ client has no configurable request timeout, so a
# stalled server-side job or network black hole can hang a query forever.
# Queries are run in a worker thread and bounded by this wall-clock timeout
# instead; a hang is treated the same as any other transient failure (see
# psr_to_gaia). The thread itself can't be killed if it's genuinely stuck --
# it's abandoned and the underlying connection will eventually time out on
# its own -- but the calling loop is freed to retry or move on.
_QUERY_EXECUTOR = ThreadPoolExecutor(max_workers=4)


def filter_position_uncertainty(df, max_arcsec=1.0):
    """Keep only pulsars with a confidently known position.

    Removes pulsars whose RA or Dec uncertainty is unmeasured (NaN) or
    exceeds ``max_arcsec``, since these can't be confidently matched against
    Gaia astrometry.

    Args:
        df (pandas.DataFrame): Pulsar table from :func:`read_atnf_long_with_errors`,
            with ``raj_err_arcsec`` and ``decj_err_arcsec`` columns.
        max_arcsec (float, optional): Maximum acceptable position uncertainty,
            in arcsec, on each axis.

    Returns:
        pandas.DataFrame: The rows that pass the cut.
    """
    mask = (
        df["raj_err_arcsec"].notna()
        & df["decj_err_arcsec"].notna()
        & (df["raj_err_arcsec"] < max_arcsec)
        & (df["decj_err_arcsec"] < max_arcsec)
    )
    return df[mask].reset_index(drop=True)


def filter_binary(df):
    """Keep only pulsars known to be in a binary.

    Args:
        df (pandas.DataFrame): Pulsar table with a ``binary_type`` column,
            where ATNF marks isolated pulsars with ``'*'``.

    Returns:
        pandas.DataFrame: The rows that pass the cut.
    """
    return df[df["binary_type"] != "*"].reset_index(drop=True)


def filter_has_proper_motion(df):
    """Keep only pulsars with a measured proper motion.

    Both propagating a pulsar's position to the Gaia epoch and confirming a
    candidate by proper-motion agreement require a measured proper motion.
    Pulsars ATNF marks as unmeasured in PMRA or PMDEC can't be used by
    either step -- attempting to propagate a NaN proper motion produces a
    NaN sky position, which Gaia's archive deterministically rejects with
    an HTTP 500 (this was found by a real full-catalogue run crashing on
    exactly this case: PSR J0045-7319, which has no measured proper motion).

    Args:
        df (pandas.DataFrame): Pulsar table with ``pmra_masyr``/``pmdec_masyr`` columns.

    Returns:
        pandas.DataFrame: The rows that pass the cut.
    """
    mask = df["pmra_masyr"].notna() & df["pmdec_masyr"].notna()
    return df[mask].reset_index(drop=True)


def _load_globular_cluster_names(gc_names_path=None):
    if gc_names_path is None:
        package_dir = os.path.dirname(os.path.abspath(__file__))
        gc_names_path = os.path.join(package_dir, "gc_pulsar_names.csv")
    with open(gc_names_path) as f:
        return {line.split()[0] for line in f if line.strip()}


def filter_in_globular(df, gc_names=None):
    """Remove pulsars known to live in a globular cluster.

    Globular clusters are crowded enough that a position/proper-motion
    cross-match against Gaia can't be trusted there.

    Args:
        df (pandas.DataFrame): Pulsar table with a ``jname`` column.
        gc_names (set, optional): Set of globular-cluster pulsar JNames to
            exclude. Defaults to the bundled ``gc_pulsar_names.csv``.

    Returns:
        pandas.DataFrame: The rows that pass the cut.
    """
    if gc_names is None:
        gc_names = _load_globular_cluster_names()
    return df[~df["jname"].isin(gc_names)].reset_index(drop=True)


def _run_cone_search(coord, radius_arcsec):
    job = Gaia.cone_search_async(coordinate=coord, radius=u.Quantity(radius_arcsec, u.arcsec))
    return job.get_results().to_pandas()


def psr_to_gaia(
    jname,
    raj_deg,
    decj_deg,
    pmra_masyr,
    pmdec_masyr,
    posepoch_mjd,
    radius_arcsec=1.0,
    max_retries=5,
    backoff_base_seconds=5.0,
    query_timeout_seconds=120.0,
):
    """Cone-searches Gaia DR3 for a possible companion to one pulsar.

    Propagates the pulsar's position from its ATNF timing epoch to the Gaia
    DR3 reference epoch (2016.0) using its proper motion, then searches for
    Gaia sources within ``radius_arcsec`` of that propagated position.

    Transient failures (network errors, an HTTP error from the Gaia archive,
    or the query hanging past ``query_timeout_seconds``) are retried with
    exponential backoff, since a large batch of queries is likely to hit at
    least one such failure.

    Args:
        jname (str): Name of the pulsar being checked for matches.
        raj_deg (float): Right ascension at ``posepoch_mjd``, in degrees.
        decj_deg (float): Declination at ``posepoch_mjd``, in degrees.
        pmra_masyr (float): Proper motion in RA, in mas/yr.
        pmdec_masyr (float): Proper motion in Dec, in mas/yr.
        posepoch_mjd (float): Epoch of ``raj_deg``/``decj_deg``, in MJD.
        radius_arcsec (float, optional): Radius of the Gaia cone search, in arcsec.
        max_retries (int, optional): Number of attempts before giving up.
        backoff_base_seconds (float, optional): Base delay for exponential
            backoff between retries (doubles each attempt).
        query_timeout_seconds (float, optional): Maximum time to wait for one
            attempt before treating it as failed and retrying. Needed because
            astroquery's Gaia client has no built-in request timeout, so a
            stalled server-side job would otherwise hang forever.

    Returns:
        pandas.DataFrame: Gaia DR3 sources found in the search (may be empty).

    Raises:
        requests.exceptions.RequestException: If every retry attempt fails.
        concurrent.futures.TimeoutError: If every retry attempt times out.
    """
    gaia_epoch = Time("2016.0", format="jyear").jyear
    p_epoch = Time(posepoch_mjd, format="mjd").jyear
    year_diff = (gaia_epoch - p_epoch) * u.yr

    pmra_degyr = (pmra_masyr * u.mas / u.yr).to(u.deg / u.yr)
    pmdec_degyr = (pmdec_masyr * u.mas / u.yr).to(u.deg / u.yr)

    new_ra = raj_deg * u.deg + pmra_degyr * year_diff
    new_dec = decj_deg * u.deg + pmdec_degyr * year_diff

    Gaia.ROW_LIMIT = 2000
    Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source"
    coord = SkyCoord(ra=new_ra, dec=new_dec, frame="icrs")

    for attempt in range(max_retries):
        try:
            future = _QUERY_EXECUTOR.submit(_run_cone_search, coord, radius_arcsec)
            return future.result(timeout=query_timeout_seconds)
        except (requests.exceptions.RequestException, FutureTimeoutError) as exc:
            if attempt == max_retries - 1:
                raise
            reason = "timed out" if isinstance(exc, FutureTimeoutError) else "failed"
            wait_seconds = backoff_base_seconds * (2**attempt)
            print(
                f"Gaia query for {jname} {reason} (attempt {attempt + 1}/{max_retries}); "
                f"retrying in {wait_seconds:.0f}s"
            )
            time.sleep(wait_seconds)


def get_matches(
    df,
    radius_arcsec=1.0,
    max_retries=5,
    backoff_base_seconds=5.0,
    query_timeout_seconds=120.0,
    checkpoint_file=None,
):
    """Cross-matches each pulsar in ``df`` against Gaia DR3.

    A pulsar whose query still fails after exhausting retries (e.g. because
    of some data problem specific to that pulsar, or a server-side hang) is
    skipped, not fatal to the whole run -- it's recorded in
    ``<checkpoint_file>.failed`` for later review, and matching continues
    with the remaining pulsars.

    Args:
        df (pandas.DataFrame): Pulsar table from :func:`read_atnf_long_with_errors`
            (optionally filtered by :func:`filter_position_uncertainty`,
            :func:`filter_binary`, :func:`filter_has_proper_motion`,
            :func:`filter_in_globular`).
        radius_arcsec (float, optional): Radius of the Gaia cone search, in arcsec.
        max_retries (int, optional): See :func:`psr_to_gaia`.
        backoff_base_seconds (float, optional): See :func:`psr_to_gaia`.
        query_timeout_seconds (float, optional): See :func:`psr_to_gaia`.
        checkpoint_file (str, optional): Path to a CSV to incrementally save
            progress to after each pulsar. If it already exists (from a
            previous run), pulsars already recorded there (successful or
            permanently failed) are skipped rather than re-queried, so an
            interrupted run can be resumed by calling this again with the
            same path.

    Returns:
        pandas.DataFrame: One row per (pulsar, Gaia candidate) pair, with all
            of the pulsar's own columns carried through alongside Gaia's.
            Empty if no pulsar had any candidates.
    """
    done_log = f"{checkpoint_file}.done" if checkpoint_file else None
    failed_log = f"{checkpoint_file}.failed" if checkpoint_file else None
    processed = set()
    frames = []

    if checkpoint_file is not None:
        if os.path.exists(done_log):
            with open(done_log) as f:
                processed = {line.strip() for line in f if line.strip()}
        if os.path.exists(checkpoint_file):
            existing = pd.read_csv(checkpoint_file)
            if len(existing) > 0:
                frames.append(existing)

    for row in df.itertuples(index=False):
        pulsar_attrs = row._asdict()
        jname = pulsar_attrs["jname"]
        if jname in processed:
            continue

        try:
            gaia_matches = psr_to_gaia(
                jname,
                pulsar_attrs["raj_deg"],
                pulsar_attrs["decj_deg"],
                pulsar_attrs["pmra_masyr"],
                pulsar_attrs["pmdec_masyr"],
                pulsar_attrs["posepoch_mjd"],
                radius_arcsec=radius_arcsec,
                max_retries=max_retries,
                backoff_base_seconds=backoff_base_seconds,
                query_timeout_seconds=query_timeout_seconds,
            )
        except (requests.exceptions.RequestException, FutureTimeoutError) as exc:
            print(f"Giving up on {jname} after {max_retries} attempts ({exc}); skipping")
            processed.add(jname)
            if checkpoint_file is not None:
                with open(done_log, "a") as f:
                    f.write(jname + "\n")
                with open(failed_log, "a") as f:
                    f.write(jname + "\n")
            continue

        if len(gaia_matches) > 0:
            for key, value in pulsar_attrs.items():
                gaia_matches[key] = value
            frames.append(gaia_matches)

        processed.add(jname)
        if checkpoint_file is not None:
            combined = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
            combined.to_csv(checkpoint_file, index=False)
            with open(done_log, "a") as f:
                f.write(jname + "\n")

    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True).drop_duplicates()


def confirm_proper_motion(matches_df, n_sigma=3.0):
    """Flags Gaia candidates whose proper motion agrees with the pulsar's.

    Position alone doesn't confirm a physical association -- a true binary
    companion should also share the pulsar's proper motion. This compares
    each candidate's Gaia proper motion against the pulsar's ATNF timing
    proper motion, combining their uncertainties in quadrature.

    Args:
        matches_df (pandas.DataFrame): Output of :func:`get_matches`.
        n_sigma (float, optional): Maximum allowed combined-uncertainty
            significance, on each axis, to count as agreement.

    Returns:
        pandas.DataFrame: ``matches_df`` with added ``pm_sigma_ra``,
            ``pm_sigma_dec``, and boolean ``pm_match`` columns.
    """
    df = matches_df.copy()

    sigma_ra = np.sqrt(df["pmra_error"] ** 2 + df["pmra_err_masyr"] ** 2)
    sigma_dec = np.sqrt(df["pmdec_error"] ** 2 + df["pmdec_err_masyr"] ** 2)

    df["pm_sigma_ra"] = np.abs(df["pmra"] - df["pmra_masyr"]) / sigma_ra
    df["pm_sigma_dec"] = np.abs(df["pmdec"] - df["pmdec_masyr"]) / sigma_dec
    df["pm_match"] = (df["pm_sigma_ra"] < n_sigma) & (df["pm_sigma_dec"] < n_sigma)
    return df


def add_gaia_distance(matches_df, min_parallax_significance=3.0):
    """Adds a Gaia-based distance estimate to each match.

    Prefers Gaia DR3's own ``distance_gspphot`` (a photometric+parallax
    distance with a proper Galactic prior) where available, falling back to
    a naive inverse-parallax distance otherwise.

    A low- or negative-significance parallax (``parallax / parallax_error``
    below ``min_parallax_significance``) makes any distance estimate for
    that source unreliable -- and is itself a sign the match may be a chance
    alignment rather than a real companion. These rows are flagged, not
    dropped: even an unreliable match is worth keeping on record to revisit
    against future Gaia data releases or other optical surveys.

    Args:
        matches_df (pandas.DataFrame): Output of :func:`get_matches` (or
            :func:`confirm_proper_motion`), with Gaia's ``parallax``,
            ``parallax_error``, and (if present) ``distance_gspphot``,
            ``distance_gspphot_lower``, ``distance_gspphot_upper``,
            ``parallax_over_error`` columns.
        min_parallax_significance (float, optional): Minimum
            ``parallax / parallax_error`` to trust a distance estimate.

    Returns:
        pandas.DataFrame: ``matches_df`` with added ``gaia_distance_pc``,
            ``gaia_distance_lower_pc``, ``gaia_distance_upper_pc``,
            ``distance_method`` (``"gspphot"`` or ``"parallax_inverse"``),
            ``parallax_significance``, and boolean
            ``low_parallax_significance`` columns.
    """
    df = matches_df.copy()
    n = len(df)

    def _column_or_nan(name):
        return df[name] if name in df.columns else pd.Series([np.nan] * n, index=df.index)

    gspphot = _column_or_nan("distance_gspphot")
    gspphot_lower = _column_or_nan("distance_gspphot_lower")
    gspphot_upper = _column_or_nan("distance_gspphot_upper")

    fallback_dist = 1000.0 / df["parallax"]
    fallback_err = np.abs(fallback_dist * df["parallax_error"] / df["parallax"])

    use_gspphot = gspphot.notna()
    df["gaia_distance_pc"] = np.where(use_gspphot, gspphot, fallback_dist)
    df["gaia_distance_lower_pc"] = np.where(use_gspphot, gspphot_lower, fallback_dist - fallback_err)
    df["gaia_distance_upper_pc"] = np.where(use_gspphot, gspphot_upper, fallback_dist + fallback_err)
    df["distance_method"] = np.where(use_gspphot, "gspphot", "parallax_inverse")

    parallax_significance = _column_or_nan("parallax_over_error").fillna(df["parallax"] / df["parallax_error"])
    df["parallax_significance"] = parallax_significance
    df["low_parallax_significance"] = parallax_significance < min_parallax_significance
    return df


def add_dm_distance(matches_df, method="ymw16"):
    """Adds a dispersion-measure-based distance estimate to each match.

    Converts each pulsar's DM to a distance using a Galactic free-electron-
    density model, for comparison against the Gaia-based distance.

    Requires the optional ``pygedm`` dependency; install it separately
    (``pip install pygedm``) if this raises ``ImportError``.

    Args:
        matches_df (pandas.DataFrame): Output of :func:`get_matches`, with
            ``gl_deg``, ``gb_deg``, and ``dm`` columns.
        method (str, optional): Electron-density model to use (``"ymw16"``
            or ``"ne2001"``).

    Returns:
        pandas.DataFrame: ``matches_df`` with an added ``dm_distance_pc`` column.
    """
    import pygedm

    df = matches_df.copy()
    distances = []
    for gl, gb, dm in zip(df["gl_deg"], df["gb_deg"], df["dm"]):
        if np.isnan(dm):
            distances.append(np.nan)
            continue
        dist, _ = pygedm.dm_to_dist(gl * u.deg, gb * u.deg, dm, method=method)
        distances.append(dist.to(u.pc).value)
    df["dm_distance_pc"] = distances
    return df


def pretty_print(match_row, filename):
    """Appends a human-readable summary of one pulsar/Gaia match to a file.

    Args:
        match_row (pandas.Series): One row from a matches DataFrame that has
            been through :func:`confirm_proper_motion` and
            :func:`add_gaia_distance` (and, optionally, :func:`add_dm_distance`).
        filename (str): Base name of the text file to append to (``.txt`` is added).
    """
    with open(filename + ".txt", "a") as g:
        g.write(f"PSR JNAME: {match_row['jname']}\n")
        g.write("\n----Gaia Match Attributes----\n")
        g.write(f"Gaia Source ID: {match_row['source_id']}\n")
        g.write(f"Gaia PMRA: {match_row['pmra']} ± {match_row['pmra_error']} mas/yr\n")
        g.write(f"Gaia PMDEC: {match_row['pmdec']} ± {match_row['pmdec_error']} mas/yr\n")
        g.write(f"Gaia Parallax: {match_row['parallax']} ± {match_row['parallax_error']} mas\n")
        g.write(f"Gaia G-band mag: {match_row['phot_g_mean_mag']}\n")
        g.write(
            f"Gaia distance ({match_row['distance_method']}): {match_row['gaia_distance_pc']:.1f} pc "
            f"[{match_row['gaia_distance_lower_pc']:.1f}, {match_row['gaia_distance_upper_pc']:.1f}]\n"
        )
        g.write("\n----Pulsar Attributes----\n")
        g.write(f"ATNF DM: {match_row['dm']} pc/cm^3\n")
        if "dm_distance_pc" in match_row:
            g.write(f"ATNF DM distance: {match_row['dm_distance_pc']:.1f} pc\n")
        g.write(f"ATNF PMRA: {match_row['pmra_masyr']} ± {match_row['pmra_err_masyr']} mas/yr\n")
        g.write(f"ATNF PMDEC: {match_row['pmdec_masyr']} ± {match_row['pmdec_err_masyr']} mas/yr\n")
        g.write(
            f"Proper motion agreement: {'YES' if match_row['pm_match'] else 'no'} "
            f"(sigma_ra={match_row['pm_sigma_ra']:.1f}, sigma_dec={match_row['pm_sigma_dec']:.1f})\n"
        )
        g.write("\n-------------------------------------------------------------\n")


def pretty_print_matches(matches_df, filename):
    """Appends a human-readable summary of every row in ``matches_df``.

    Args:
        matches_df (pandas.DataFrame): See :func:`pretty_print`.
        filename (str): Base name of the text file to append to (``.txt`` is added).
    """
    for _, row in matches_df.iterrows():
        pretty_print(row, filename)


def matching_pipeline(
    input_file,
    output_file,
    max_pos_err_arcsec=1.0,
    radius_arcsec=1.0,
    n_sigma=3.0,
    min_parallax_significance=3.0,
    gc_names=None,
    include_dm_distance=True,
    pretty_print_output=False,
    max_retries=5,
    backoff_base_seconds=5.0,
    query_timeout_seconds=120.0,
    checkpoint_file=None,
):
    """Runs the full ATNF-to-Gaia cross-match pipeline end to end.

    Reads an ATNF "long with errors" export, filters to pulsars with a
    confident position that are in a binary and not in a globular cluster,
    cross-matches the survivors against Gaia DR3, confirms candidates by
    proper-motion agreement, and adds Gaia- and DM-based distance estimates.

    Args:
        input_file (str): Path to the ATNF "long with errors" text export.
        output_file (str): Path to write the final matches CSV to.
        max_pos_err_arcsec (float, optional): See :func:`filter_position_uncertainty`.
        radius_arcsec (float, optional): See :func:`get_matches`.
        n_sigma (float, optional): See :func:`confirm_proper_motion`.
        min_parallax_significance (float, optional): See :func:`add_gaia_distance`.
        gc_names (set, optional): See :func:`filter_in_globular`.
        include_dm_distance (bool, optional): Whether to add a DM-based
            distance via :func:`add_dm_distance` (requires ``pygedm``).
        pretty_print_output (bool, optional): Whether to also write a
            human-readable ``output_file + '.txt'`` via :func:`pretty_print_matches`.
        max_retries (int, optional): See :func:`get_matches`.
        backoff_base_seconds (float, optional): See :func:`get_matches`.
        query_timeout_seconds (float, optional): See :func:`get_matches`.
        checkpoint_file (str, optional): See :func:`get_matches` -- recommended
            for any large run, so it can be resumed if interrupted.

    Returns:
        pandas.DataFrame: The final matches table (also written to ``output_file``).
    """
    df = read_atnf_long_with_errors(input_file)
    df = filter_position_uncertainty(df, max_arcsec=max_pos_err_arcsec)
    df = filter_binary(df)
    df = filter_has_proper_motion(df)
    df = filter_in_globular(df, gc_names=gc_names)

    matches = get_matches(
        df,
        radius_arcsec=radius_arcsec,
        max_retries=max_retries,
        backoff_base_seconds=backoff_base_seconds,
        query_timeout_seconds=query_timeout_seconds,
        checkpoint_file=checkpoint_file,
    )
    if len(matches) > 0:
        matches = confirm_proper_motion(matches, n_sigma=n_sigma)
        matches = add_gaia_distance(matches, min_parallax_significance=min_parallax_significance)
        if include_dm_distance:
            matches = add_dm_distance(matches)

    matches.to_csv(output_file, index=False)

    if pretty_print_output and len(matches) > 0:
        pretty_print_matches(matches, output_file)

    return matches
