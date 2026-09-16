# PSRmatch
<img alt="GitHub release (latest by date including pre-releases)" src="https://img.shields.io/github/v/release/AnnikaDeutsch/Binary-Pulsar-Distances?display_name=tag&include_prereleases"> <img alt="GitHub" src="https://img.shields.io/github/license/AnnikaDeutsch/Binary-Pulsar-Distances"> <img alt="codeastro" src="https://img.shields.io/badge/June%202022-codeastro-blueviolet">

<!-- ABOUT THE PROJECT -->
## About PSRmatch
This package, PSRmatch, cross matches pulsars with Gaia to identify pulsar binary companions that are typically
detectable in the Gaia (optical) wavelengths. It does this by taking in a list of pulsars, running them through
a number of criteria to look at only those potentially useful in pulsar timing (position uncertainty, binary
status, not in a globular cluster), searching in a small radius around the positions of each of the remaining
pulsars (propagated to the Gaia epoch using each pulsar's proper motion), and gives back a list of any Gaia DR3
sources found in that search. These sources are potentially binary companions to some of the original input
pulsars.

The end goal is to use confirmed companions to improve distance estimates for the pulsars, since Gaia parallaxes
are often better constrained than pulsar-timing or dispersion-measure-based distances. **Project status:**
the full pipeline is implemented: position-based cross-matching against Gaia DR3, confirming candidates by
proper-motion agreement (`confirm_proper_motion`), and comparing Gaia-based distances
(`add_gaia_distance`) against DM-based distances (`add_dm_distance`, via the `pygedm` YMW16 model). Producing
and writing up results across a full pulsar sample is in progress. The longer-term goal is to generalize the
pipeline to cross-match other pulsar catalogues (e.g. the MeerKAT Thousand Pulsar Array) against other optical
surveys (e.g. PanSTARRS, OGLE), not just ATNF against Gaia. See `CLAUDE.md` for the full current status and
roadmap.

<p align="right">(<a href="#readme-top">back to top</a>)</p>

### Installation

1. Clone the repo
   ```sh
   git clone https://github.com/AnnikaDeutsch/Binary-Pulsar-Distances.git
   cd Binary-Pulsar-Distances
   ```
2. Install dependencies
   ```sh
   pip install -r requirements.txt
   ```
3. Install the package (editable, so local changes are picked up immediately)
   ```sh
   pip install -e .
   ```

**Note on `pygedm`:** it compiles a C extension (the YMW16 Galactic
electron-density model) and can be finicky. If `pip install pygedm` fails to
build, or imports but raises `ModuleNotFoundError: No module named
'pkg_resources'` or `AttributeError: module 'scipy.integrate' has no
attribute 'simps'`, try installing it in a dedicated environment with
`pip install "setuptools<81" "scipy<1.14" pygedm` (newer `setuptools`/`scipy`
removed APIs `pygedm` still relies on). `add_dm_distance()` is the only
function in this package that needs `pygedm`; everything else works without it.

<p align="right">(<a href="#readme-top">back to top</a>)</p>

## Documentation

Documentation is built with Sphinx and is not currently hosted anywhere public. To build and view it locally:
```sh
pip install sphinx sphinx-rtd-theme
cd docs
make html
open _build/html/index.html
```

<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE.txt` for more information.

<p align="right">(<a href="#readme-top">back to top</a>)</p>

