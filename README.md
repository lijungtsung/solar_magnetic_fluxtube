# Solar magnetic fluxtube

This repository contains numerical methods and results for calculating solar-surface magnetic flux structures undergoing magneto-hydrostatic equilibrium with the ambient gas. The numerical methods implemented follow the approach described by Steiner, Pneuman, and Stenflo (1986).

## Contents
- **'run1_network.ipynb'**: Provides numerical results for finite-difference calculations on a 65×65 mesh grid.

- **'run2_Br.py'** and **'run2_Bz.py'**: Perform linear interpolation of the radial (Br) and axial (Bz) components of the magnetic field between grid points.

- **'run3_B_continuum.ipynb'**: Demonstrates the use of a continuum magnetic field function.

## Methodology

The numerical methods used in this project are based on:

- Steiner, O., Pneuman, G., Stenflo, J., "Numerical models for solar magnetic fluxtubes", <a href="https://ui.adsabs.harvard.edu/abs/1986A%26A...170..126S/abstract"> Astronomy and Astrophysics, 170, 126 (1986)

## Citation

If you use this code or any part of this repository in your research, please cite the following papers:

- Steiner, O., Pneuman, G. W., Stenflo, J. O., "Numerical models for solar magnetic fluxtubes", <a href="https://ui.adsabs.harvard.edu/abs/1986A%26A...170..126S/abstract"> Astronomy and Astrophysics 170, 126 (1986)

- Li, J-T, Beacom, J. F., Griffith, S., Peter, A. H. G., "Small-Scale Magnetic Fields are Critical to Shaping Solar Gamma-Ray Emission", <a href="https://doi.org/10.3847/1538-4357/ad158f"> The Astrophysical Journal 961, 167 (2024)
