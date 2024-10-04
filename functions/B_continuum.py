'''
This scripts should only be executed after "run2_Bfields.ipynb" has been executed.
'''

import numpy as np
from parameters import input_values as iv
from functions import B_lattice_final as B_LT_final


## read B_lattice data ##
N_iter_final = np.load(iv.B_LT_output_dir+'N_iter.npy')
tube_grid_point_final = np.load(iv.B_LT_output_dir+'tube_grid_point/tube_grid_point_%i.npy'%(N_iter_final) )

#############################################
################ B continuum ################
#############################################
init_BGrid = np.loadtxt(iv.B_CTM_output_dir+'BGrid_input_values.txt', delimiter='\t')
delta_z, delta_r, N_grid_z, N_grid_r = init_BGrid

Br_output = np.load(iv.B_CTM_output_dir+'BrGrid_read.npy').astype(np.float32)
Bz_output = np.load(iv.B_CTM_output_dir+'BzGrid_read.npy').astype(np.float32)

#############################################
################ tube radius ################
#############################################
heigh_global = np.arange(0, len(tube_grid_point_final), 1) * iv.hz + iv.z_bottom    ## [cm]
radius_xy = tube_grid_point_final * iv.hr                                           ## [cm]

heigh_global = np.concatenate(([-30*1e8], heigh_global, [20*1e8]), axis=0)          ## [cm]
radius_xy = np.concatenate(([radius_xy[0]], radius_xy, [radius_xy[-1]]), axis=0)    ## [cm]

if np.all(np.diff(heigh_global) >= 0) != True:
    print('Warning: heigh_global is not monotonically increasing array')

'''
"func_tube_radius" gives tube radius as a function of global z
z=0 is surface of the Sun
'''
## from scipy import interpolate
## func_tube_radius = interpolate.interp1d(heigh_global, radius_xy)
def func_tube_radius(z_arr_input):
    '''
    Numpy interpolation is faster than scipy interp1d;
    Note that numpy interp only accept monotonically increasing sample points.
    'heigh_global' must be a monotonically increasing array.
    '''
    return np.interp(z_arr_input, heigh_global, radius_xy)  ## return unit is [cm]


def B_vec(r_vec_arr):
    
    x, y, z_global = r_vec_arr
    r = np.sqrt(x**2 + y**2)
    
    if z_global < iv.z_bottom:
        z_global = iv.z_bottom 
    
    if z_global > iv.z_top:
        z_global = iv.z_top
    
    z_tube = z_global - iv.z_bottom

    j_float = round(z_tube/delta_z)
    i_float = round(r/delta_r)
    
    if i_float >= N_grid_r:
        return np.zeros(3)
    
    Bz = Bz_output[j_float][i_float]
    
    if r == 0.0:
        return np.array([0, 0, Bz])
    
    Br_over_r = Br_output[j_float][i_float] / r
    Bx = Br_over_r * x
    By = Br_over_r * y
    
    return np.array([Bx, By, Bz])   ## return unit is [Gauss]



def B_vec_interp(r_vec_arr):
    
    x, y, z_global = r_vec_arr
    r = np.sqrt(x**2 + y**2)
    
    if z_global < iv.z_bottom:
        z_global = iv.z_bottom 
    
    if z_global > iv.z_top:
        z_global = iv.z_top
        
    if r >= N_grid_r * delta_r:
        return np.zeros(3)
    
    Bz = B_LT_final.Bz_continuous(r, z_global)
    
    if r == 0.0:
        return np.array([0, 0, Bz])
    
    Br_over_r = B_LT_final.Br_continuous(r, z_global) / r
    Bx = Br_over_r * x
    By = Br_over_r * y
    
    return np.array([Bx, By, Bz])   ## return unit is [Gauss]