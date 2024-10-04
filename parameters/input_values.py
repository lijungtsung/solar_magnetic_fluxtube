'''
Input values for network field from z=600 km to z=10 Mm.
'''

import numpy as np

###
B_LT_output_dir = 'data1_output/'   ## B_lacttice output data
B_CTM_output_dir = 'data2_output/'  ## B_continuum output data

###
'''
Size of the fluxtube: two dimension -- r and z
"z_bottom" and "z_top" are global coordinates with z=0 at the solar surface
(i.e., z=0 is where tau_5000 = 1)
"L" is the vertical length of the simulation box (i.e., from bottom to top)
"W" is the horizontal length of the simulation box
'''
z_bottom = 0.7 * 1e8      # [cm], z = 700 km
z_top    = 10 * 1e8       # [cm], z = 10 Mm
W = 7 * 1e8               # [cm], W = 10 Mm
L = (z_top - z_bottom)    # [cm], L = 9.3 Mm

###
'''
The minimum amd maximum heights for GCR trajectory in the network field simulation
'''
HeightMin = z_bottom           # [cm]
HeightMax = z_top + 0.5e8      # [cm]

###
'''
Grid points
Equal number of grid points in x and y axes
'''
N_side = 2**6 + 1         # N_side is the number of grid points on each of the axis
N_spacing = N_side - 1
hz = L/N_spacing
hr = W/N_spacing
Nx_matrix = N_side-2      # -2 because left and right edges follow the Dirichlet boundary condition
Ny_matrix = N_side-1      # -1 because bottom edge follows the Dirichlet boundary condition

###
'''
Fluxtube at the solar surface
'''
B_0z = 130                     ## [gauss], Bz at the bottom of the simulation box (i.e., Bz at z_bottom)
R_star = 3 * 1e8               ## [cm], R_base = 3 Mm at z=700 km