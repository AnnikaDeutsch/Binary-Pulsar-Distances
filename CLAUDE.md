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

1. **Finish the original pipeline**: add the proper-motion-agreement
   confirmation step, then use confirmed matches to actually produce improved
   distance estimates and write up results (the `galactic_projections.ipynb`
   plotting work is a first step in this direction).
2. **Generalize the framework**: rather than being hardcoded to ATNF + Gaia,
   support cross-matching *any* pulsar catalogue (adding the **MeerKAT
   Thousand Pulsar Array**) against *any* optical survey (adding
   **PanSTARRS** and **OGLE** alongside Gaia).

## Current state of the code (updated after Phase 0 stabilization, 2026-09-16)

- `Binary_Pulsar_Distances/eliminating.py` is the maintained pipeline module:
  `check_pos_uncertainty` → `check_binary` → `check_in_globular` →
  `get_matches` (Gaia DR3 cone search per pulsar) → `matching_pipeline`
  (wraps all of the above). No proper-motion-agreement confirmation step
  exists yet — matching stops at "any Gaia source within the search radius."
- Earlier prototypes/duplicates of the position-match idea (`match_gaia_to_psr.py`
  + notebook, `oldquery.ipynb`, `query.ipynb`, and the old `matching_test.py`
  suite that tested them) have been moved to `legacy/` — see `legacy/README.md`.
  They predate `eliminating.py` and are not run or maintained, but the
  `psrqpy`-based live-ATNF-query idea in `match_gaia_to_psr.py` may be worth
  reviving when generalizing to other catalogues.
- `Binary_Pulsar_Distances/gc_names.py` + `gc_pulsars.csv` +
  `gc_pulsar_names.csv` build the globular-cluster exclusion list used by
  `check_in_globular`.
- `galactic_projections.ipynb` (+ `temp.csv`/`tempdr3.csv`) is the newest
  results work: plotting galactic positions of DR2 vs DR3 match results.
- `eliminating_test.py` now imports the real `Binary_Pulsar_Distances.eliminating`
  module (no more inline duplication) and uses paths relative to the repo,
  so it runs on any machine. Current status: `check_pos_uncertainty` and
  `check_in_globular` tests pass; `check_binary`'s test is `xfail` (real bug,
  see below); the full `matching_pipeline` test is `skip`ped (too slow to run
  routinely — see below). `get_matches`/`psr_to_gaia` have no unit test
  coverage yet (the old tests for them, in `legacy/matching_test.py`, targeted
  a superseded function signature).
- `docs/` builds cleanly with Sphinx (`cd docs && make html`) and its
  autodoc page for `Binary_Pulsar_Distances.eliminating` is in sync with the
  code. `docs/get_matches.rst` (an orphaned, broken page for a module that
  never existed standalone) was removed.
- `Binary_Pulsar_Distances/__init__.py` now exports the public API
  (`check_binary`, `check_pos_uncertainty`, `check_in_globular`, `get_matches`,
  `matching_pipeline`, `pretty_print`, `psr_to_gaia`).
- `setup.py` previously had `packages=find_packages()` commented out, so
  `pip install .` installed *no* code at all. Fixed to actually package
  `Binary_Pulsar_Distances` and pull `install_requires` from
  `requirements.txt`; verified with a clean-venv install.
- README installation/doc instructions are now accurate to the current setup
  (`pip install -e .`, local Sphinx build instructions instead of the old
  WSL-only doc link).

## Known open issues (found during Phase 0, not yet fixed)

- **`check_binary()` has a real column-indexing bug**, not just a stale test
  fixture: it reads the binary-status flag at `values[11]`, which only lines
  up with the column layout that `check_pos_uncertainty()` *outputs*
  (index, jname, ra, dec, pmra, pmra_err, pmdec, pmdec_err, posepoch, DM,
  DM_err, binary, trailing). The plain ATNF-export layout used by the
  `check_binary` unit test fixture (`text_files_test/e1_input.csv`) only has
  9 real fields with the binary flag at index 8, so calling `check_binary()`
  on it raises `IndexError`. Whether `check_pos_uncertainty`'s output format
  is actually the intended real-world ATNF "long with errors" layout (i.e.
  `check_binary` is correct and the test fixture is wrong), or `check_binary`
  should instead work on the plain layout, needs Annika's domain knowledge of
  the true ATNF export format before fixing — don't just guess indices.
  Tracked as an `xfail` in `eliminating_test.py::TestBinaryCheck`.
- The full `matching_pipeline()` integration test
  (`eliminating_test.py::TestMatchingPipeline`) is `skip`ped: run against the
  full `all_atnf.csv` fixture (~4000 pulsars), it does a live Gaia cone
  search per surviving pulsar and did not finish in over 2 minutes. Needs a
  small hand-picked fixture (a handful of pulsars) to be a routine test.
- No tests exist yet for `get_matches()`/`psr_to_gaia()` (the Gaia-querying
  functions) against the current function signatures.

## Known rough edges to keep in mind

- Several near-duplicate implementations of the same position-matching logic
  used to exist across files; the superseded ones now live in `legacy/`, but
  keep new work in `eliminating.py` (or its Phase 2 successor modules) rather
  than adding another parallel variant.
- No pluggable module/interface layout yet for the "generalize to any
  catalogue/survey" goal — `eliminating.py` is still ATNF+Gaia-specific.

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
