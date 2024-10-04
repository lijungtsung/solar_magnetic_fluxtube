import numpy as np
from parameters import input_values as iv


def Br_lattice(i, j, n_iteration):

    if n_iteration <= 0:
        print('n_iteration need to be greater or equal to 1')
        return 'n_iteration need to be greater or equal to 1'
    
    if type(i)==float or type(j)==float:
        print('i, j must be integer')
        return 'i, j must be integer'
    
    if j >= iv.N_side or j <=- 1 or i >= iv.N_side or i <=- 1:
        print('index out of bound')
        return 'index out of bound'
    
    if i == 0 or i == iv.N_side-1:
        return 0.0
        
    if j == iv.N_side-1:
        ## Neumann boundary condition. At the top, B is z direction, so Br=0.
        return 0.0
    
    # if j == 0 and i >= int(iv.R_star/iv.hr)+1:
    #     return 0.0
    
    if j >= 1 and j <= iv.N_side-2 and i >= 1 and i <= iv.N_side-2:
        r = i * iv.hr
        idx_ij = iv.N_side*i + j
        idx_up = idx_ij + 1
        idx_down = idx_ij - 1
        data = np.load(iv.B_LT_output_dir+'data/data%i.npy'%(n_iteration-1))
        return -1./r * (data[idx_up][2] - data[idx_down][2]) / (2.0*iv.hz)  # central difference
    
    if j == 0 and i >= 1 and i <= iv.N_side-2:
        r = i * iv.hr
        idx_ij = iv.N_side*i + j
        idx_up = idx_ij + 1
        data = np.load(iv.B_LT_output_dir+'data/data%i.npy'%(n_iteration-1))
        return -1./r * (data[idx_up][2] - data[idx_ij][2]) / (1.0*iv.hz)    # forward difference 


############################################################        
def Bz_lattice(i, j, n_iteration):

    if n_iteration <= 0:
        print('n_iteration need to be greater or equal to 1')
        return 'n_iteration need to be greater or equal to 1'
    
    if type(i)==float or type(j)==float:
        print('i, j must be integer')
        return 'i, j must be integer'
    
    if j >= iv.N_side or j <= -1 or i >= iv.N_side or i <= -1:
        print('index out of bound')
        return 'index out of bound'
    
    if i == iv.N_side-1:
        return 0.0
    
    if i == 0:
        # Apply forward difference
        r = 0.5 * iv.hr
        idx_ij = iv.N_side*i + j
        idx_right = idx_ij + iv.N_side     # N_side*(i+1) + j
        data = np.load(iv.B_LT_output_dir+'data/data%i.npy'%(n_iteration-1))
        return 1./r * (data[idx_right][2] - data[idx_ij][2]) / (1.0*iv.hr)       # forward difference
    
    if i <= iv.N_side-2 and i >= 1:
        # Apply central difference
        r = i * iv.hr
        idx_ij = iv.N_side*i + j
        idx_right = idx_ij + iv.N_side     # N_side*(i+1) + j
        idx_left  = idx_ij - iv.N_side     # N_side*(i-1) + j
        data = np.load(iv.B_LT_output_dir+'data/data%i.npy'%(n_iteration-1))
        return 1./r * (data[idx_right][2] - data[idx_left][2]) / (2.0*iv.hr)     # central difference

############################################################     
def B_lattice(i, j, n_iteration):
    if n_iteration <= 0:
        print('n_iteration need to be greater or equal to 1')
        return 'n_iteration need to be greater or equal to 1'
    
    if n_iteration >= 1:
        return np.sqrt(Br_lattice(i, j, n_iteration)**2 + Bz_lattice(i, j, n_iteration)**2)    
    
    
def Br_external_lattice(i, j, n_iteration):

    if n_iteration <= 0:
        print('n_iteration need to be greater or equal to 1')
        return 'n_iteration need to be greater or equal to 1'
    
    if type(i)==float or type(j)==float:
        print('i, j must be integer')
        return 'i, j must be integer'
    
    if j >= iv.N_side or j<=-1 or i >= iv.N_side or i<=-1:
        print('index out of bound')
        return 'index out of bound'
    
    if i == 0 or i == iv.N_side-1:
        return 0.0
        
    if j == iv.N_side-1:
        return 0.0
    
    if j == 0 and i >= 1 and i <= iv.N_side-2:
        
        tube_grid_point = np.load(iv.B_LT_output_dir+'tube_grid_point/tube_grid_point_%i.npy'%(n_iteration) )
        if i == tube_grid_point[j]:
            i = i+1
        return Br_lattice(i, j, n_iteration)  
        
    if j >= 1 and j <= iv.N_side-2 and i >= 1 and i <= iv.N_side-2:
        
        tube_grid_point = np.load(iv.B_LT_output_dir+'tube_grid_point/tube_grid_point_%i.npy'%(n_iteration) )
        if i == tube_grid_point[j]:
            i = i+1
        
        LR_list = np.load(iv.B_LT_output_dir+'nearby_point/Lower_Right%i.npy'%(n_iteration)).tolist()
        if [i, j] in LR_list:
            # Apply backward difference
            r = i * iv.hr
            idx_ij   = iv.N_side*i + j
            idx_down = idx_ij   - 1
            data = np.load(iv.B_LT_output_dir+'data/data%i.npy'%(n_iteration-1))
            return -1./r * (data[idx_ij][2] - data[idx_down][2]) / (1.0*iv.hz)    # backward difference
        else:
            return Br_lattice(i, j, n_iteration)
    
    
############################################################        
def Bz_external_lattice(i, j, n_iteration):

    if n_iteration <= 0:
        print('n_iteration need to be greater or equal to 1')
        return 'n_iteration need to be greater or equal to 1'
    
    if type(i)==float or type(j)==float:
        print('i, j must be integer')
        return 'i, j must be integer'
    
    if j >= iv.N_side or j <= -1 or i >= iv.N_side or i<=-1:
        print('index out of bound')
        return 'index out of bound'
    
    if i == iv.N_side-1:
        return 0.0   
    
    if i <= iv.N_side-2 and i >= 0:
        
        tube_grid_point = np.load(iv.B_LT_output_dir+'tube_grid_point/tube_grid_point_%i.npy'%(n_iteration) )
        if i == tube_grid_point[j]:
            i = i+1
        
        r = (i+0.5) * iv.hr 
        idx_ij    = iv.N_side*i + j
        idx_right = idx_ij   + iv.N_side
        data = np.load(iv.B_LT_output_dir+'data/data%i.npy'%(n_iteration-1))
        return 1./r * (data[idx_right][2] - data[idx_ij][2]) / (1.0*iv.hr)    # forward difference

############################################################     
def B_external_lattice(i, j, n_iteration):
    if n_iteration <= 0:
        print('n_iteration need to be greater or equal to 1')
        return 'n_iteration need to be greater or equal to 1'
    
    if n_iteration >= 1:
        return np.sqrt(Br_external_lattice(i, j, n_iteration)**2 + Bz_external_lattice(i, j, n_iteration)**2)    
    
def Br_internal_lattice(i, j, n_iteration):

    if n_iteration <= 0:
        print('n_iteration need to be greater or equal to 1')
        return 'n_iteration need to be greater or equal to 1'
    
    if type(i)==float or type(j)==float:
        print('i, j must be integer')
        return 'i, j must be integer'
    
    if j >= iv.N_side or j <= -1 or i >= iv.N_side or i <= -1:
        print('index out of bound')
        return 'index out of bound'
    
    if i == 0 or i == iv.N_side-1:
        return 0.0
        
    if j == iv.N_side-1:
        return 0.0
    
    if j == 0 and i >= 1 and i <= iv.N_side-2:
        
        tube_grid_point = np.load(iv.B_LT_output_dir+'tube_grid_point/tube_grid_point_%i.npy'%(n_iteration) )
        if i == tube_grid_point[j]:
            i = i-1   
        return Br_lattice(i, j, n_iteration)  
        
    if j >= 1 and j <= iv.N_side-2 and i >= 1 and i <= iv.N_side-2:
        
        tube_grid_point = np.load( iv.B_LT_output_dir+'tube_grid_point/tube_grid_point_%i.npy'%(n_iteration) )
        if i == tube_grid_point[j]:
            i = i-1        
        
        UL_list = np.load(iv.B_LT_output_dir+'nearby_point/Upper_Left%i.npy'%(n_iteration)).tolist()
        if [i, j] in UL_list:
            r = i * iv.hr
            idx_ij = iv.N_side*i + j
            idx_up = idx_ij   + 1
            data = np.load(iv.B_LT_output_dir+'data/data%i.npy'%(n_iteration-1))
            return -1./r * (data[idx_up][2] - data[idx_ij][2]) / (1.0*iv.hz)    # forward difference
        else:
            return Br_lattice(i, j, n_iteration)

############################################################        
def Bz_internal_lattice(i, j, n_iteration):

    if n_iteration <= 0:
        print('n_iteration need to be greater or equal to 1')
        return 'n_iteration need to be greater or equal to 1'
    
    if type(i)==float or type(j)==float:
        print('i, j must be integer')
        return 'i, j must be integer'
    
    if j >= iv.N_side or j<=-1 or i >= iv.N_side or i<=-1:
        print('index out of bound')
        return 'index out of bound'
    
    if i == iv.N_side-1:
        return 0.0
    
    if i <= iv.N_side-2 and i >= 0:
        
        tube_grid_point = np.load( iv.B_LT_output_dir+'tube_grid_point/tube_grid_point_%i.npy'%(n_iteration) )
        if i == tube_grid_point[j]:
            i = i-1
        
        r = (i-0.5) * iv.hr
        idx_ij       = iv.N_side*i + j
        idx_backward = idx_ij   - iv.N_side
        data = np.load(iv.B_LT_output_dir+'data/data%i.npy'%(n_iteration-1))
        return 1./r * (data[idx_ij][2] - data[idx_backward][2]) / (1.0*iv.hr)     # Backward difference
    
############################################################     
def B_internal_lattice(i, j, n_iteration):
    if n_iteration <= 0:
        print('n_iteration need to be greater or equal to 1')
        return 'n_iteration need to be greater or equal to 1'
    
    if n_iteration >= 1:
        return np.sqrt(Br_internal_lattice(i, j, n_iteration)**2 + Bz_internal_lattice(i, j, n_iteration)**2)      
############################################################    