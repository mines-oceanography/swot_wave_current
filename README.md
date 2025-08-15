
# Source code for 
Villas Bôas, A. B., Marechal, G., and Bohé. Observing Interactions Between Waves, Winds, and Currents from SWOT. Geophysical Research Letters, Accepted.

# Abstract
The Surface Water and Ocean Topography (SWOT) satellite mission enables, for the first
time, two-dimensional (2D) mapping of significant wave height ($H_s$) at kilometer-scale
resolution. Using data from SWOT’s Ka-band Radar Interferometer (KaRIn), this study
investigates interactions between surface waves, winds, and currents across diverse dynamic regimes, including western boundary currents, mesoscale turbulence, tropical cyclones, and wave group modulation. SWOT reveals unprecedented 2D spatial gradients
in $H_s$, capturing fine-scale variability previously identified only in numerical models. These
observations highlight the critical role of currents in shaping the wave field and show strong
agreement with theoretical predictions. SWOT’s high-resolution wave data represent a
transformative advance in understanding $H_s$ gradients, paving the way for refining operational models and addressing challenges in characterizing the influence of sea state
gradients on coupled air-sea processes.

# Authors
* [Bia Villas Boas](https://mines-oceanography.github.io/) <villasboas@mines.edu>
* [Gwendal Marechal](https://gmarechal.github.io/)
* [Alejandro Bohé](https://ieeexplore.ieee.org/author/37089945952)

# Data
All processed data and model output used in this manuscript are available in the `data` directory in this repository. 

# Funding
This project was funded by NASA through the [SWOT](https://swot.jpl.nasa.gov/) program (award 80NSSC24K1640), the Ocean Vector Winds Science Team (80NSSC23K0979), and the Earth System Explorer program (award 80GSFC24CA067) supporting development of the [ODYSEA](https://odysea.ucsd.edu/) Concept Study Report. Additional support was obtained from the [ONR MURI program](https://www.minesnewsroom.com/news/mines-researcher-flying-eye-storm-learn-more-about-air-sea-interactions) (grant N00014-24-1-2554). 

# How to use this repository
All figures in the manuscript can be reproduced using the Python scripts from this repository and the [data](https://github.com/mines-oceanography/swot_wave_current/tree/main/data). To do so, follow these steps

1. Make a local copy of this repository by either cloning or downloading it.

Your directory tree should look like this:

```
swot_wave_current/
├── data/
├── figures/
├── notebooks/
├── src/
└── pixi.toml
```
2. If you would like to reproduce the figures from the manuscript using this repository, we recommend that you install [Pixi](https://pixi.sh/latest/). You can then run:

```
$ pixi shell
```
from the project root. This will build the dependencies and activate the Python environment for this project. You can then launch Jupyter Lab as you would normally:

```
$ jupyter lab
```

If you follow the steps above, you should be able to reproduce all figures by running the notebooks from the notebooks directory without having to adjust any paths.

# How to cite this code
DOI coming soon.
