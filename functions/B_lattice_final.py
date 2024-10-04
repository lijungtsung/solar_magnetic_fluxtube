import numpy as np
from scipy import interpolate
from parameters import input_values as iv
from functions import B_lattice as B_LT


N_iter_final = np.load(iv.B_LT_output_dir+'N_iter.npy')
data_final = np.load(iv.B_LT_output_dir+'data/data%i.npy'%(N_iter_final))
tube_grid_point_final = np.load(iv.B_LT_output_dir+'tube_grid_point/tube_grid_point_%i.npy'%(N_iter_final) )
UL_list_final = np.load(iv.B_LT_output_dir+'nearby_point/Upper_Left%i.npy'%(N_iter_final)).tolist()
LR_list_final = np.load(iv.B_LT_output_dir+'nearby_point/Lower_Right%i.npy'%(N_iter_final)).tolist()


def Bz_final_lattice(i, j, n_iteration):
    
    if n_iteration <= 0:
        print('n_iteration need to be greater or equal to 1')
        return 'n_iteration need to be greater or equal to 1'
    
    if i > int(tube_grid_point_final[j])+1:
        return B_LT.Bz_lattice(i, j, n_iteration)
    
    if i == int(tube_grid_point_final[j])+1:
        return B_LT.Bz_external_lattice(i, j, n_iteration)
    
    if i == tube_grid_point_final[j]:
        return B_LT.Bz_internal_lattice(i-1, j, n_iteration)
    
    if i == int(tube_grid_point_final[j]):
        return B_LT.Bz_internal_lattice(i, j, n_iteration)
    
    if i < int(tube_grid_point_final[j]):
        return B_LT.Bz_lattice(i, j, n_iteration)
    
def Br_final_lattice(i, j, n_iteration):
    
    if n_iteration <= 0:
        print('n_iteration need to be greater or equal to 1')
        return 'n_iteration need to be greater or equal to 1'
    
    if i > int(tube_grid_point_final[j])+1:
        return B_LT.Br_lattice(i, j, n_iteration)
    
    if i == int(tube_grid_point_final[j])+1:
        return B_LT.Br_external_lattice(i, j, n_iteration)
    
    if i == tube_grid_point_final[j]:
        return B_LT.Br_internal_lattice(i-1, j, n_iteration)
    
    if i == int(tube_grid_point_final[j]):
        return B_LT.Br_internal_lattice(i, j, n_iteration)
    
    if i < int(tube_grid_point_final[j]):
        return B_LT.Br_lattice(i, j, n_iteration)
    
def B_final_lattice(i, j, n_iteration):
    
    if n_iteration <= 0:
        print('n_iteration need to be greater or equal to 1')
        return 'n_iteration need to be greater or equal to 1'
    
    return np.sqrt(Bz_final_lattice(i, j, n_iteration)**2 + Br_final_lattice(i, j, n_iteration)**2)



############################################################
###############   B final, continuous   ####################
############################################################


def Bz_continuous(r, z_global):

    if r >= iv.W:
        return 0.
    
    z_tube = z_global - iv.z_bottom
    i_float, j_float = r/iv.hr, z_tube/iv.hz
    
    if z_global >= iv.z_top:
        B00 = Bz_final_lattice(int(i_float), int(j_float), N_iter_final)
        B01 = Bz_final_lattice(int(i_float)+1, int(j_float), N_iter_final)
        return np.interp(i_float, np.array([int(i_float), int(i_float)+1]),  np.array([B00, B01]))
    
    i_arr = np.array([int(i_float), int(i_float)+1])
    j_arr = np.array([int(j_float), int(j_float)+1])
    M00 = Bz_final_lattice(i_arr[0], j_arr[0], N_iter_final)
    M01 = Bz_final_lattice(i_arr[1], j_arr[0], N_iter_final)
    M10 = Bz_final_lattice(i_arr[0], j_arr[1], N_iter_final)
    M11 = Bz_final_lattice(i_arr[1], j_arr[1], N_iter_final)
    matrix = np.array([[M00, M01], [M10, M11]])
    intp_2d = interpolate.interp2d(i_arr, j_arr, matrix, kind='linear')
    
    return intp_2d(i_float, j_float)[0]


def Br_continuous(r, z_global):
    
    if r >= iv.W:
        return 0.
    
    z_tube = z_global - iv.z_bottom
    i_float, j_float = r/iv.hr, z_tube/iv.hz
    
    if z_global >= iv.z_top:
        B00 = Br_final_lattice(int(i_float), int(j_float), N_iter_final)
        B01 = Br_final_lattice(int(i_float)+1, int(j_float), N_iter_final)
        return np.interp(i_float, np.array([int(i_float), int(i_float)+1]),  np.array([B00, B01]))
    
    i_arr = np.array([int(i_float), int(i_float)+1])
    j_arr = np.array([int(j_float), int(j_float)+1])
    M00 = Br_final_lattice(i_arr[0], j_arr[0], N_iter_final)
    M01 = Br_final_lattice(i_arr[1], j_arr[0], N_iter_final)
    M10 = Br_final_lattice(i_arr[0], j_arr[1], N_iter_final)
    M11 = Br_final_lattice(i_arr[1], j_arr[1], N_iter_final)
    matrix = np.array([[M00, M01], [M10, M11]])
    intp_2d = interpolate.interp2d(i_arr, j_arr, matrix, kind='linear')
    
    return intp_2d(i_float, j_float)[0]