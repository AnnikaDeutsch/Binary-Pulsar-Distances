# PSRmatch
<img alt="GitHub release (latest by date including pre-releases)" src="https://img.shields.io/github/v/release/AnnikaDeutsch/Binary-Pulsar-Distances?display_name=tag&include_prereleases"> <img alt="GitHub" src="https://img.shields.io/github/license/AnnikaDeutsch/Binary-Pulsar-Distances"> <img alt="codeastro" src="https://img.shields.io/badge/June%202022-codeastro-blueviolet">

<!-- ABOUT THE PROJECT -->
## About PSRmatch

This project began when I participated in a one week long workshop called codeastro. In this workshop, we
had the objective of building a small software package. The package I began in that week is the project I 
continued to work on for the rest of the summer.  

This package, PSRmatch, cross matches pulsars with Gaia to identify pulsar binary companions that are typically
detectable in the Gaia (optical) wavelengths. It does this by taking in a list of pulsars, running them through
a number of criteria to look at only those potentially useful in pulsar timing, searching in a small radius
around the positions of each of the remaining pulsars, and gives back a list of any Gaia sources to result
from this search. These sources are potentially binary companions to some of the original input pulsars. 
Further checks can be performed on these output matches, comparing the astrometric parameters of the Gaia sources
to known parameters of certain pulsar companions, such as G-band magnitude. 

<p align="right">(<a href="#readme-top">back to top</a>)</p>

### Installation

1. Clone the repo
   ```sh
   git clone https://github.com/AnnikaDeutsch/Binary-Pulsar-Distances.git
   ```
3. Install dependencies
   ```sh
   pip install -r requirements.txt
   ```
4. Install the package
   ```js
   pip install psrmatch
   ```

<p align="right">(<a href="#readme-top">back to top</a>)</p>

## Documentation

Documentation can be accessed here: file://wsl.localhost/Ubuntu/home/annika_deutsch/Binary-Pulsar-Distances/docs/_build/html/index.html

<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE.txt` for more information.

<p align="right">(<a href="#readme-top">back to top</a>)</p>

