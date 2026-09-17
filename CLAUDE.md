# PSRmatch (Binary-Pulsar-Distances)

## NEVER MAKE COMMITS

**Do not run `git commit` (or `git push`) in this repo under any
circumstances, even if asked to "commit" as part of a larger task or if it
seems like the natural next step after finishing work.** Annika (the user)
always makes commits herself. Stage or leave changes as working-tree edits
and let her review and commit them.

## Origin and goal

This project started as Annika Deutsch's undergraduate research project. The
core scientific idea: cross-match the **ATNF pulsar catalogue**
(https://www.atnf.csiro.au/research/pulsar/psrcat/) against **Gaia DR3**
(https://gea.esac.esa.int/archive/) to find optical companions to binary
pulsars. Matching is done first by **sky position** (propagated to the Gaia
epoch using each pulsar's proper motion) and should then be confirmed by
**proper-motion agreement** between the pulsar and the candidate Gaia source.
A confirmed match gives an optical counterpart with a Gaia parallax/distance,
which is generally better constrained than pulsar-timing or DM-based distance
estimates — the payoff is improved distance estimates for the pulsars in
question.

The project was left unfinished at the position-matching stage. It is now
being resumed with two goals, in order:

1. **Finish the original pipeline** (done — see below): add the
   proper-motion-agreement confirmation step, then use confirmed matches to
   actually produce improved distance estimates. Producing and writing up
   results across a full pulsar sample (the `galactic_projections.ipynb`
   plotting work is a first step in this direction) is still in progress.
2. **Generalize the framework**: rather than being hardcoded to ATNF + Gaia,
   support cross-matching *any* pulsar catalogue (adding the **MeerKAT
   Thousand Pulsar Array**) against *any* optical survey (adding
   **PanSTARRS** and **OGLE** alongside Gaia).

## Current state of the code (updated after Phase 1, 2026-09-16)

- The pipeline was rearchitected around a pandas DataFrame with named,
  unit-normalized columns instead of raw `;`/whitespace-delimited rows
  addressed by positional index. This was necessary, not just a cleanup:
  proper-motion agreement needs ATNF's PM uncertainties to survive from the
  raw catalogue export all the way to the Gaia-comparison step, and the old
  row format truncated them away partway through. It also resolved the
  Phase-0-deferred `check_binary` indexing bug by construction (no more
  competing index conventions once columns are named).
- `Binary_Pulsar_Distances/atnf.py` (new): `read_atnf_long_with_errors()`
  parses ATNF psrcat's "long with errors" export format into that DataFrame,
  including a real unit conversion for RAJ_ERR (ATNF gives it in seconds of
  time; converted to arcsec via `* 15 * cos(dec)` before comparing against
  DECJ_ERR, which ATNF already gives in arcsec) and Galactic `gl_deg`/`gb_deg`
  (derived via `astropy.coordinates.SkyCoord`, needed for DM-distance).
- `Binary_Pulsar_Distances/eliminating.py` is now: `filter_position_uncertainty`
  (real arcsec cutoff, replacing the old exponent-sign heuristic) →
  `filter_binary` → `filter_has_proper_motion` (see below) →
  `filter_in_globular` → `get_matches` (Gaia DR3 cone search per pulsar,
  carrying every ATNF column through onto each matched Gaia row) →
  `confirm_proper_motion` (flags `pm_match` via a 3-sigma
  combined-uncertainty comparison of ATNF vs. Gaia proper motion) →
  `add_gaia_distance` (prefers Gaia DR3's own `distance_gspphot`, falls back
  to naive inverse-parallax; flags `low_parallax_significance` -- see below)
  → `add_dm_distance` (DM → distance via the `pygedm` YMW16 model).
  `matching_pipeline()` orchestrates all of the above. The old
  `check_binary`/`check_pos_uncertainty`/`check_in_globular`/raw
  `get_matches` functions were replaced outright, not kept as shims.
- `filter_has_proper_motion()`: **found via a real full-catalogue run**
  crashing on PSR J0045-7319, which has no measured ATNF proper motion at
  all. Propagating a NaN PM produces a NaN sky position, which Gaia's
  archive deterministically rejects with an HTTP 500 every time -- not
  transient, so retries just wasted time before still crashing. Turned out
  **116 of the 255 pulsars surviving the other filters** (nearly half) lack
  a measured PM, so this filter now runs before `get_matches()`
  (255 → 139 candidates actually queried). A pulsar without PM can't be
  propagated *or* PM-confirmed anyway, so this isn't a big loss of
  information -- it's information the pipeline literally can't use yet.
- `add_gaia_distance()`'s `low_parallax_significance` flag: **found in the
  first real full run's output** -- PSR J1543-5149's match had
  `parallax_over_error ≈ -0.06` (consistent with noise, not a real
  parallax), which the naive inverse-parallax fallback turned into a
  physically meaningless **-56089 pc**. Per Annika's call, these rows are
  **flagged, not dropped** (`low_parallax_significance` boolean column,
  threshold `min_parallax_significance=3.0` on
  `parallax_over_error`/computed `parallax/parallax_error`) -- kept on
  record to revisit against future Gaia data releases or other optical
  surveys, since a low-significance match now isn't necessarily wrong,
  just unconfirmable with current data.
- Removed a hardcoded special case in the old `psr_to_gaia` that forced a
  fixed 5 arcsec search radius for one specific pulsar (`J0437-4715`) with no
  documented justification — inconsistent with the "generalizable" goal, and
  callers can already pass their own `radius_arcsec`.
- `pygedm` **is verified working**, in a dedicated conda env named `pygedm`
  (`/opt/anaconda3/envs/pygedm`, currently Python 3.10). Run anything needing
  `pygedm` with `/opt/anaconda3/envs/pygedm/bin/python -m pytest
  eliminating_test.py` (that env also has the rest of `requirements.txt`).
  **This took real effort to get working and is fragile — read this before
  touching that env again:**
  - The root problem is that **this machine's Xcode Command Line Tools
    currently point at a broken/malformed SDK** (`MacOSX27.0.sdk`, referenced
    under `/Library/Developer/CommandLineTools/SDKs/`), which breaks *any*
    fresh compile of `pygedm`'s C/C++ extensions (YMW16, NE2001), regardless
    of Python version — confirmed by watching a from-source build fail
    identically on both Python 3.14 and a freshly recreated Python 3.10 env.
    This is **not** a pygedm-version or Python-version issue; earlier
    guessing that "just use Python 3.10" would fix it was wrong. The only
    reason anything works at all is that pip has a **locally cached
    prebuilt wheel** from an earlier successful compile (in
    `~/Library/Caches/pip/wheels/...pygedm-3.3.0-cp310-cp310-macosx_11_0_arm64.whl`
    — this cache is keyed to Python 3.10 specifically and is **user-level,
    not tied to any one conda env**, so it survives `conda env remove`).
    Until the system SDK is fixed, **don't run `pip install --no-binary
    pygedm` or otherwise force a rebuild** — it will fail. Stick to Python
    3.10 in this env so pip keeps reusing that cached wheel.
  - Fresh `setuptools` (>=81) drops `pkg_resources`, which `pygedm/__init__.py`
    still imports → `pip install "setuptools<81"` in that env.
  - Fresh `scipy` (>=1.14) renamed `integrate.simps` to `integrate.simpson`,
    and `pygedm`'s YMW16 wrapper still calls the old name → `pip install
    "scipy<1.14"` in that env.
  - The cached wheel's NE2001 extension (`ne21c`) links against
    `libf2c.dylib`, which isn't available anywhere on this machine (not in
    conda-forge or Homebrew) — **patched around it** by editing the
    installed `.../site-packages/pygedm/pygedm.py` to wrap `from . import
    ne2001_wrapper` in a `try/except ImportError` (falls back to `None`).
    We only ever use `method='ymw16'`, so this is fine; calling
    `dm_to_dist(..., method='ne2001')` in that env will now fail instead
    (not something we do). **This is a hand-patch to installed
    site-packages, not tracked by git or pip — it will be silently lost if
    `pygedm` is ever reinstalled/upgraded in that env**, and would need to
    be reapplied (same edit) if that happens.
  - Confirmed `pygedm.dm_to_dist(gl, gb, dm, method='ymw16')` matches my
    assumed API exactly, and `add_dm_distance()` reproduces the published
    VLBI distance to PSR J2222-0137 (267.3 pc, Guo et al. 2021) to within
    0.1 pc using its real ATNF DM — verified twice, independently, across
    the original and recreated envs.
  - The real, permanent fix would be repairing this machine's Xcode Command
    Line Tools/SDK (out of scope to do unilaterally — it's a system-wide
    change, not project-scoped); until then, treat the `pygedm` env as a
    fragile, hand-assembled artifact rather than something safely
    reproducible by just re-running `pip install -r requirements.txt`.
- Earlier prototypes/duplicates of the position-match idea (`match_gaia_to_psr.py`
  + notebook, `oldquery.ipynb`, `query.ipynb`, and the old `matching_test.py`
  suite that tested them) live in `legacy/` — see `legacy/README.md`. The
  `psrqpy`-based live-ATNF-query idea in `match_gaia_to_psr.py` may be worth
  reviving when generalizing to other catalogues; `psrqpy` was also used
  ad hoc (not as a runtime dependency) to fetch real parameters for
  `text_files_test/small_atnf_known_binary.csv` (see tests, below).
- `Binary_Pulsar_Distances/gc_names.py` + `gc_pulsars.csv` +
  `gc_pulsar_names.csv` build the globular-cluster exclusion list used by
  `filter_in_globular`.
- `galactic_projections.ipynb` (+ `temp.csv`/`tempdr3.csv`) predates the
  Phase 1 rearchitecture and plots the *old* pipeline's DR2/DR3 match output
  format — it will need updating once we produce results from the new
  pipeline (Phase 1's "produce and write up results" item, still open).
- `eliminating_test.py` was rewritten against the new DataFrame-based API:
  unit tests for `read_atnf_long_with_errors`, `filter_position_uncertainty`,
  `filter_binary`, `filter_in_globular`, `confirm_proper_motion`, and
  `add_gaia_distance` all pass without network access. `TestMatchingPipeline`
  now runs the *full* pipeline end-to-end (live Gaia query) against a tiny,
  real fixture (one real pulsar, J2222-0137, a well-known confirmed binary —
  see `text_files_test/small_atnf_known_binary.csv`) instead of being
  permanently skipped; it's a plumbing/smoke test, not a check that this
  specific companion is Gaia-detectable (many white-dwarf companions are too
  faint for Gaia, so `pm_match` may legitimately be `False` for every
  candidate — confirmed this is in fact the case for the 5 Gaia sources
  found near J2222-0137 within 60 arcsec). `TestAddDmDistance` checks
  `add_dm_distance()` against the published VLBI distance for the same
  pulsar; it `pytest.importorskip("pygedm")`s, so it's skipped wherever
  `pygedm` isn't installed and runs for real in the `pygedm` conda env.
  `TestPsrToGaiaRetry` and `TestGetMatchesCheckpointing` cover the new
  retry/checkpoint behavior via `monkeypatch` (no real network flakiness or
  interrupted runs needed to test it).
- `docs/` builds cleanly with Sphinx (`cd docs && make html`); added
  `docs/atnf.rst` alongside `docs/eliminating.rst`, both wired into
  `docs/index.rst`'s toctree.
- `Binary_Pulsar_Distances/__init__.py` exports the current public API
  (`read_atnf_long_with_errors`, `filter_position_uncertainty`,
  `filter_binary`, `filter_in_globular`, `get_matches`,
  `confirm_proper_motion`, `add_gaia_distance`, `add_dm_distance`,
  `matching_pipeline`, `pretty_print`, `pretty_print_matches`, `psr_to_gaia`).
- **Runtime cost of a full run, measured empirically (not estimated):** a
  live Gaia `cone_search_async` query averages **~21.5s** per pulsar (timed
  over 6 real queries; one of the 8 attempted hit a live `HTTP 500` from the
  Gaia archive). Running the real filter chain against `all_atnf.csv`
  (3320 valid rows) leaves **255 pulsars** to actually query (1844 pass
  position, 292 of those are binaries, 255 of those aren't in a globular
  cluster) — so a full run is roughly 255 × ~21.5s ≈ **~90 minutes**, serial,
  scaling with the *filtered ATNF pulsar count*, not with Gaia DR3's size
  (each query is a cheap indexed positional lookup against Gaia's archive,
  regardless of its ~1.8 billion sources).
- `psr_to_gaia()`/`get_matches()` now have **retry-with-exponential-backoff**
  (`max_retries`, `backoff_base_seconds`, catching
  `requests.exceptions.RequestException`) and `get_matches()` supports a
  `checkpoint_file` that incrementally saves progress and a `.done` sibling
  log of completed JNames, so an interrupted run (e.g. a persistent Gaia
  outage after retries are exhausted) can be resumed by calling it again
  with the same `checkpoint_file` path rather than restarting from scratch.
  A pulsar that still fails after exhausting retries is skipped (not fatal
  to the whole batch) and recorded in `<checkpoint_file>.failed` for review.
  `matching_pipeline()` passes all of this through.
- **`psr_to_gaia()` also has a `query_timeout_seconds` (default 120s) hard
  wall-clock timeout per attempt, found necessary by a real run hanging.**
  `astroquery`'s Gaia TAP+ client has no configurable request timeout at
  all -- a real run sat with an `ESTABLISHED` TCP connection to
  `gea.esac.esa.int` and 0% CPU for **over an hour** with zero progress,
  because the underlying HTTP call simply never returned or errored (so the
  existing retry-on-exception logic never even triggered). Each attempt is
  now run in a worker thread (`concurrent.futures.ThreadPoolExecutor`,
  module-level `_QUERY_EXECUTOR`) bounded by `future.result(timeout=...)`;
  a timeout is treated as just another retryable failure. Note this can't
  actually kill a truly stuck request -- the abandoned thread lingers until
  its own connection eventually resolves -- but it stops the *caller* from
  blocking forever. (In that same real run, the stuck query happened to
  resolve on its own after roughly an hour, and the run finished
  successfully before this fix was even deployed -- but don't count on that;
  the timeout is the actual fix.)
- **Investigated whether a single batched query (upload the whole pulsar
  list to Gaia's TAP+ service, one ADQL join against `gaiadr3.gaia_source`)
  would be faster than one query per pulsar — empirically, it was not.**
  Tested with the same 8 pulsars used for the per-query timing above: the
  batch upload-crossmatch job took **862.9 seconds** (~14 minutes) for all
  8, versus ~172s doing them one at a time. The fixed per-query overhead
  clearly isn't what dominates; the join itself seems to be the expensive
  part on Gaia's end (possibly costed as a full join with no index pushed
  down onto the uploaded table, or a lower-priority job queue for
  user-table-upload jobs — not confirmed which). **Don't switch to a batch
  upload-based approach without re-testing** — this could change if
  Gaia's archive/query planner behaves differently at a larger n, but there's
  no evidence for that yet, so per-pulsar querying (now with retry/checkpoint)
  remains the recommended approach. The CDS X-Match service
  (`astroquery.xmatch`) was considered but not tested as an alternative.

## First full real run (2026-09-16/17)

Ran `matching_pipeline()` against the full `all_atnf.csv` fixture (via the
`pygedm` env, `checkpoint_file='full_crossmatch_checkpoint.csv'`). Output:
`full_crossmatch_results.csv` (final, with `pm_match`/distances/significance
flag) and `full_crossmatch_checkpoint.csv` (raw pre-distance matches, so
distance/significance logic can be recomputed later without re-querying
Gaia -- this is how the `low_parallax_significance` flag was retroactively
added to this run's output).

- **139/139 pulsars processed, 0 permanent failures.** Took much longer than
  the small-sample estimate (see the ~90-minute estimate above) — closer to
  a few hours, mostly due to the hour-long hang described above plus a
  slower-than-sampled average query time; the ~21.5s/query figure was from a
  small sample and shouldn't be taken as a tight bound.
- **24 pulsars** have at least one Gaia candidate within 1 arcsec; **13 pass
  proper-motion confirmation** (`pm_match=True`).
- **Of those 13, 7 have `low_parallax_significance=True`** (J0337+1715,
  J1048+2339, J1543-5149, J1653-0158, J1959+2048, J2039-5617, J2339-0533) —
  so only **6** (J0437-4715, J1012+5307, J1023+0038, J1227-4853, J1417-4402,
  J1723-2837) currently have both PM agreement *and* a trustworthy distance.
  This is the real headline number for "how many pulsars got a usable
  improved distance from this run," not the 13.
- **Sanity checks against known sources both passed:** PSR J0437-4715 (a
  famous nearby MSP; the pulsar itself is the historical Gaia detection,
  not a companion) came back at 140.8 pc via inverse-parallax vs. its own
  YMW16 DM-distance of 156.1 pc — good agreement, and also confirms
  removing the old hardcoded 5-arcsec special case for this exact pulsar
  (see above) didn't break anything. PSR J1023+0038 (a well-known
  transitional MSP) came back at 1417.9 pc via `gspphot` vs. 1112.4 pc via
  DM — same ballpark, independent methods.
- galactic_projections.ipynb has not yet been updated to plot this output
  (still targets the old pre-Phase-1 format) -- next real step.

## Known open issues

- `pygedm` still doesn't compile in the main working environment on this
  machine (see above) — only the dedicated `pygedm` conda env has it. Not an
  issue on Annika's side; just remember to use that env for anything
  DM-distance-related until/unless the main env's toolchain is fixed.
- `galactic_projections.ipynb` targets the pre-Phase-1 output format/columns
  and needs updating before it can plot new pipeline output (now that a real
  results file, `full_crossmatch_results.csv`, actually exists to plot).
- A `pd.concat` `FutureWarning` shows up in `get_matches()`'s checkpoint
  writes (empty/all-NA frame concatenation) -- cosmetic, not yet cleaned up.

## Known rough edges to keep in mind

- Several near-duplicate implementations of the same position-matching logic
  used to exist across files; the superseded ones now live in `legacy/`, but
  keep new work in `eliminating.py`/`atnf.py` (or their Phase 2 successor
  modules) rather than adding another parallel variant.
- No pluggable module/interface layout yet for the "generalize to any
  catalogue/survey" goal — `atnf.py`/`eliminating.py` are still ATNF+Gaia-
  specific, though the DataFrame-based design from Phase 1 should make that
  refactor more mechanical (a catalogue reader just needs to populate the
  same column schema; a survey backend just needs to consume it and return
  Gaia-shaped candidate columns).

## Working agreement: keep documentation current

**Whenever code changes in this repo (new/renamed functions, new modules,
changed usage/CLI, new catalogue or survey backends), update in the same
pass:**

- `docs/*.rst` — add or update autodoc stubs so every public module/function
  is documented, and fix `docs/index.rst`'s toctree accordingly.
- `README.md` — keep the "About," installation, and usage sections accurate
  to the current package layout and capabilities.

Don't defer this to a follow-up — treat doc/README updates as part of
finishing any code change here, not an optional extra.
