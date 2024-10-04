# Solar magnetic fluxtube

In this repository, I show the numerical method for calculating the magnetic flux structures undergoing magneto-hydrostatic equilibrium with the ambient gas. The numerical method I follow can be found in Steiner, Pneuman, and Stenflo 1986 (A&A 170, 126) <https://ui.adsabs.harvard.edu/abs/1986A%26A...170..126S/abstract>

run1_network.ipynb provides the numerical result for finite-difference calculation in a 65 * 65 mesh grid. 

run2_Br.py and run2_Bz.py calculate the linear interpolation of Br and Bz between the grid points.

run3_B_continuum.ipynb provides the demonstration for using continuum magnetic field function.
