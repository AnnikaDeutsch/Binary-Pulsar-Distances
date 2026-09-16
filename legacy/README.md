# Legacy / superseded code

These files are early prototypes from the original undergrad project that
predate `Binary_Pulsar_Distances/eliminating.py`, which is the maintained
pipeline implementation. They're kept for reference (some ideas here, like
querying ATNF live via `psrqpy` instead of a pre-exported text file, may be
worth reviving) but are not run, tested, or maintained.

- `oldquery.ipynb`, `query.ipynb` — near-duplicate exploratory notebooks;
  the original prototyping ground for the position-propagation + Gaia-query
  approach that became `eliminating.py`.
- `match_gaia_to_psr.py`, `match_gaia_to_psr.ipynb` — an alternate
  implementation that queries ATNF live via `psrqpy` (requires the `psrqpy`
  package, not in the active `requirements.txt`) and does a Gaia box search
  (`Gaia.query_object_async`) instead of a cone search, keyed on pulsar name
  or an explicit RA/Dec/radius.
- `matching_test.py` — a pytest suite for an older `psr_to_gaia`/`get_matches`
  signature (`height`/`width` params, Gaia DR2) that no longer matches
  `eliminating.py`'s current functions, and references `.vot.gz` fixture
  files that aren't present in the repo. Superseded by the tests in
  `eliminating_test.py`.
