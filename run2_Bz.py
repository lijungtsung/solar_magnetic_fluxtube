import numpy as np
import os
from parameters import run2_parameters as param
from parameters import input_values as iv
from functions import B_lattice_final as B_LT_final
from mpi4py.futures import MPIPoolExecutor


delta_z = iv.hz/param.Nz_intp    # [cm]
delta_r = iv.hr/param.Nr_intp    # [cm]
N_grid_z = int(iv.L/delta_z) + 1
N_grid_r = int(iv.W/delta_r) + 1
    

f = open(iv.B_CTM_output_dir+'BGrid_input_values.txt', 'w')
f.write('# delta_z = hz/%i, [cm]\n'%param.Nz_intp)
f.write('# delta_r = hr/%i, [cm]\n'%param.Nr_intp)
f.write('# N_grid_z = int(L/delta_z) + 1\n')
f.write('# N_grid_r = int(W/delta_r) + 1\n')
f.write('# \n')
f.write('# Columns: delta_z[cm], delta_r[cm], N_grid_z, N_grid_r\n')
f.write('%e\t%e\t%i\t%i'%(delta_z, delta_r, N_grid_z, N_grid_r))
f.close()




def Bz_ij_grid(i_float, j_float):
    z_tube = j_float * delta_z
    z_here_global = z_tube + iv.z_bottom
    r_here = i_float * delta_r
    return B_LT_final.Bz_continuous(r_here, z_here_global)


def wrapping_calculate_Bz_row_z(idx_z):
    
    '''
    i_float = idx_r, j_float = idx_z.
    '''
    Bz_row_z_arr = []
    for idx_r in range(N_grid_r):
        Bz_row_z_arr += [Bz_ij_grid(idx_r, idx_z)]
        
    Bz_row_z_arr = np.array(Bz_row_z_arr, dtype=np.float32)
    np.save(iv.B_CTM_output_dir+'temp_Bz_rows/row%i.npy'%idx_z, Bz_row_z_arr)    

    return None


#################################################
###########  run parallel computing  ############
#################################################

#################################################
if __name__ == '__main__':
    
    
    #### check whether temp_Bz_rows folders exist or not ####
    directory_path = os.path.join(os.getcwd(), iv.B_CTM_output_dir, 'temp_Bz_rows')
    os.makedirs(directory_path, exist_ok=True)
    del directory_path
    
    
    #### start the multi-processing ####
    with MPIPoolExecutor() as pool:
        pool.map(wrapping_calculate_Bz_row_z, range(N_grid_z))
    
    
    #### join all the rows ####
    Bz_all_rows = [None] * N_grid_z
    
    for idx_z in range(N_grid_z):
        temp_arr = np.load(iv.B_CTM_output_dir+'temp_Bz_rows/row%i.npy'%idx_z)
        Bz_all_rows[idx_z] = temp_arr.tolist()

    Bz_all_rows = np.array(Bz_all_rows, dtype=np.float32)
    np.save(iv.B_CTM_output_dir+'BzGrid_read.npy', Bz_all_rows)
    del Bz_all_rows
    
    
    #### delete files in the 'temp_rows' ####
    directory_path = os.path.join(os.getcwd(), iv.B_CTM_output_dir, 'temp_Bz_rows')
    files = os.listdir(directory_path)
    for file in files:
        file_path = os.path.join(directory_path, file)
        if os.path.isfile(file_path):
            os.remove(file_path)    
    #### end ####