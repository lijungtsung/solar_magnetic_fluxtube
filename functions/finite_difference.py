import numpy as np

from functions import PhotosphereModel as PM
from parameters import input_values as iv
from functions import B_lattice as B_LT


############################
#### boundary condition ####
############################

def coef_R(i): # i+1
    return iv.L/iv.W - iv.L/(2.0*i*iv.W)
def coef_L(i): # i-1
    return iv.L/iv.W + iv.L/(2.0*i*iv.W)
def coef_T(j): # j+1
    return iv.W/iv.L
def coef_B(j): # j-1
    return iv.W/iv.L

coef_ctr = -2.0 * (iv.W/iv.L + iv.L/iv.W)

def psi_star(i):
    r = i * iv.hr
    if r <= iv.R_star:
        return 0.5 * iv.B_0z * r**2
    if r > iv.R_star:
        return 0.5 * iv.B_0z * iv.R_star**2
    
## left edge
boundary_0j = np.zeros(iv.N_side) # left, j = [0, j_last]
## bottom edge
boundary_i0 = np.array([psi_star(i) for i in range(iv.N_side)])
## right edge
boundary_last_j = np.array([psi_star(iv.N_side)]*iv.N_side)  # right, j = [0, j_last]  
####################################################################################



#####################################
#### calculate matrix and source ####
#####################################

def calculate_matrix(n_iteration):
    '''
    This function calculates the matrix elements for the elliptical equation (the Grad-Shafranov equation).
    The size of the matrix is (Nx_matrix*Ny_matrix, Nx_matrix*Ny_matrix)
    '''
    
    ## creating empty matrix and souce vector
    matrix = np.zeros(( iv.Nx_matrix*iv.Ny_matrix, iv.Nx_matrix*iv.Ny_matrix ))
    np.fill_diagonal( matrix, coef_ctr )
    source = calculate_source_cs(n_iteration)
    
    #### Calculate scalar function u at each grid point. 
    #### Scalar function u at the grid points at the boundaries are not calculated here.
    for i in range(1, iv.Nx_matrix+1):
        for j in range(1, iv.Ny_matrix):
            
            ## the index of (i,j) for the source array
            source_idx = (i-1)*iv.Ny_matrix + j - 1   # checked

            if i==1:
                source[source_idx] += -1 * coef_L(i) * boundary_0j[j]
                matrix_idx = (i+1-1)*iv.Ny_matrix + j - 1
                matrix[source_idx][matrix_idx] = coef_R(i)

            elif i==iv.Nx_matrix:
                source[source_idx] += -1 * coef_R(i) * boundary_last_j[j]
                matrix_idx = (i-1-1)*iv.Ny_matrix + j - 1
                matrix[source_idx][matrix_idx] = coef_L(i)

            else:
                matrix_idx = (i-1-1)*iv.Ny_matrix + j - 1
                matrix[source_idx][matrix_idx] = coef_L(i)
                matrix_idx = (i+1-1)*iv.Ny_matrix + j - 1
                matrix[source_idx][matrix_idx] = coef_R(i)
            
            if j==1:
                source[source_idx] += -1 * coef_B(j) * boundary_i0[i]
                matrix_idx = (i-1)*iv.Ny_matrix + (j+1) - 1
                matrix[source_idx][matrix_idx] = coef_T(j)

            else:
                matrix_idx = (i-1)*iv.Ny_matrix + (j-1) - 1
                matrix[source_idx][matrix_idx] = coef_B(j)
                matrix_idx = (i-1)*iv.Ny_matrix + (j+1) - 1
                matrix[source_idx][matrix_idx] = coef_T(j)

    #### Calculate scalar function u at the top of the boundary (Neumann boundary condition)
    j_last = iv.Ny_matrix
    for i in range(1, iv.Nx_matrix+1):
        source_idx = (i-1) * iv.Ny_matrix + j_last - 1 

        if i==1:
            source[source_idx] += -1 * coef_L(i) * boundary_0j[j_last]
            matrix_idx = (i+1-1)*iv.Ny_matrix + j_last - 1
            matrix[source_idx][matrix_idx] = coef_R(i)
            matrix_idx = (i-1)*iv.Ny_matrix + (j_last-1) - 1
            matrix[source_idx][matrix_idx] = coef_B(j_last) + coef_T(j_last)

        elif i==iv.Nx_matrix:
            source[source_idx] += -1 * coef_R(i) * boundary_last_j[j_last]
            matrix_idx = (i-1-1)*iv.Ny_matrix + j_last - 1
            matrix[source_idx][matrix_idx] = coef_L(i)
            matrix_idx = (i-1)*iv.Ny_matrix + (j_last-1) - 1
            matrix[source_idx][matrix_idx] = coef_B(j_last) + coef_T(j_last)

        else:
            matrix_idx = (i-1-1)*iv.Ny_matrix + j_last - 1
            matrix[source_idx][matrix_idx] = coef_L(i)
            matrix_idx = (i+1-1)*iv.Ny_matrix + j_last - 1
            matrix[source_idx][matrix_idx] = coef_R(i)
            matrix_idx = (i-1)*iv.Ny_matrix + (j_last-1) - 1
            matrix[source_idx][matrix_idx] = coef_B(j_last) + coef_T(j_last)
    
    return matrix, source


def calculate_u_sol(n_iteration, theta_relaxation):
    
    matrix_n, source_n = calculate_matrix(n_iteration)
    
    if n_iteration == 0:
        u_sol_n = np.linalg.inv(matrix_n).dot(source_n)
        np.save(iv.B_LT_output_dir+'u_sol/u_sol%i.npy'%(n_iteration), u_sol_n)
        return u_sol_n
    
    if n_iteration >= 1:
        '''
        Use the implicit under-relaxation method
        For under ralaxation, theta_relaxation < 1.0
        '''
        for i in range(len(matrix_n)):
            matrix_n[i][i] = matrix_n[i][i]/theta_relaxation
        
        u_sol_old = np.load(iv.B_LT_output_dir+'u_sol/u_sol%i.npy'%(n_iteration-1))
        source_n = source_n + (1.0-theta_relaxation)/theta_relaxation * coef_ctr * u_sol_old
        u_sol_n = np.linalg.inv(matrix_n).dot(source_n)
        
        np.save(iv.B_LT_output_dir+'u_sol/u_sol%i.npy'%(n_iteration), u_sol_n)
        return u_sol_n

def calculate_data(n_iteration, theta_relaxation):
    '''
    This function assembles the final solution for u. 
    The mesh points at the boundaries r=0, r=W, and z_bottom=0 are included.
    data[0] corresponds to (i=0, j=0) with data structure as [hr*0, hz*0, boundary_0j[0]],
    data[1] corresponds to (i=0, j=1) with data structure as [hr*0, hz*1, boundary_0j[1]],
    data[2] corresponds to (i=0, j=2) with data structure as [hr*0, hz*2, boundary_0j[1]],
    etc...
    '''
    u_sol = calculate_u_sol(n_iteration, theta_relaxation)
    data = []
    
    for i in range(iv.N_side):
        for j in range(iv.N_side):
            
            ## at the left boundary
            if i==0:
                data += [[iv.hr*i, iv.hz*j, boundary_0j[j]]]
            
            ## at the right boundary
            if i==iv.N_side-1:
                data += [[iv.hr*i, iv.hz*j, boundary_last_j[j]]]

            ## the bottom boundary, but exclude the most left and right mesh points
            if i!=0 and i!=iv.N_side-1 and j==0:
                data += [[iv.hr*i, iv.hz*j, boundary_i0[i]]]
            
            ## the mesh points not at the boundaries
            if i!=0 and i!=iv.N_side-1 and j!=0:
                source_idx = (i-1) * iv.Ny_matrix + j - 1 
                data += [[iv.hr*i, iv.hz*j, u_sol[source_idx][0]]]   
    
    np.save(iv.B_LT_output_dir+'data/data%i.npy'%n_iteration, data)  
    return data

def Delta_c(n_iteration):
    u_sol_n = np.load(iv.B_LT_output_dir+'u_sol/u_sol%i.npy'%(n_iteration))
    u_sol_n_previous = np.load(iv.B_LT_output_dir+'u_sol/u_sol%i.npy'%(n_iteration-1))
    difference_sum = np.sum( ((u_sol_n - u_sol_n_previous)/u_sol_n_previous)**2 )
    return np.sqrt( difference_sum / (iv.Nx_matrix*iv.Ny_matrix) )
####################################################################################



###############################
#### calculate tube radius ####
###############################

def phi_radius(j, n_iteration, flux):    
    '''
    This function calculates the radius of the tube with flux = flux.
    The input flux should be smaller than the total flux of the entire fluxtube.
    '''
    
    ## first, calculate total vertical fluxes at the bottom of the tube
    total_magnetic_flux = np.pi * (0.5*iv.hr)**2 * iv.B_0z   # magnetic flux at i=0
    for i in range(1, int(iv.R_star/iv.hr)):
        total_magnetic_flux += 2 * np.pi * i * iv.hr**2 * iv.B_0z 
    total_magnetic_flux += np.pi * iv.hr**2 * iv.B_0z * 0.5 *  ( 2*int(iv.R_star/iv.hr) - 0.5 )
    
    ## Next, calculate flux at the height z; j = z/hz
    off_i = 0.001
    magnetic_flux_count = np.pi * (0.5*iv.hr)**2 * B_LT.Bz_lattice(0, j, n_iteration)
    for i_tube_radius in range(1, iv.N_side-1):
        d_flux_i = 2 * np.pi * i_tube_radius * iv.hr**2 * B_LT.Bz_lattice(i_tube_radius, j, n_iteration)
        magnetic_flux_count += d_flux_i
        
        if d_flux_i == 0:
            print('Error: d_flux_i = 0')
            return 'error'
        
        if magnetic_flux_count >= flux:
            delta_flux = flux - (magnetic_flux_count - d_flux_i)
            i_tube_edge = np.sqrt( (i_tube_radius-0.5)**2  
                        + delta_flux/(np.pi*iv.hr**2 * B_LT.Bz_internal_lattice(i_tube_radius, j, n_iteration)) )
            if i_tube_edge >= iv.N_side-2:
                return iv.N_side-2-off_i
            else: 
                return i_tube_edge
        
        if magnetic_flux_count < flux and i_tube_radius == (iv.N_side-2):
            return iv.N_side-2-off_i

################
def tube_radius(j, n_iteration):
    '''
    This function calculates the radius of the flux tube. 
    The output is R_tube/hr where R_tube is the radius of the tube at height j*hz from the bottom of tube.
    '''
    
    if n_iteration == 0:
        return 'n_iteration should be a postive integer'
    
    ## first, calculate the total flux at the base 
    magnetic_flux = np.pi * (0.5*iv.hr)**2 * iv.B_0z 
    for i in range(1, int(iv.R_star/iv.hr)):
        magnetic_flux += 2*np.pi * i * iv.hr**2 * iv.B_0z 
    magnetic_flux += np.pi * iv.hr**2 * iv.B_0z * 0.5 *  ( 2*int(iv.R_star/iv.hr) - 0.5 )
    
    ## Next, calculate flux at the height z; j = z/hz
    off_i = 0.001
    magnetic_flux_count = np.pi * (0.5*iv.hr)**2 * B_LT.Bz_lattice(0, j, n_iteration)
    for i_tube_radius in range(1, iv.N_side-1):
        d_flux_i = 2.0 * np.pi * i_tube_radius * iv.hr**2 * B_LT.Bz_lattice(i_tube_radius, j, n_iteration)
        magnetic_flux_count += d_flux_i
        
        if d_flux_i == 0:
            return 'error'
        
        elif magnetic_flux_count >= magnetic_flux:
            delta_flux = magnetic_flux - (magnetic_flux_count - d_flux_i) 
            i_tube_edge = np.sqrt( (i_tube_radius-0.5)**2  
                        + delta_flux/(np.pi*iv.hr**2*B_LT.Bz_lattice(i_tube_radius, j, n_iteration) ) )
            if i_tube_edge >= iv.N_side-2:
                return iv.N_side-2-off_i
            else: 
                return i_tube_edge
        
        elif magnetic_flux_count < magnetic_flux and i_tube_radius == (iv.N_side-2):
            return iv.N_side-2-off_i
        
        else:
            None    ## continue the loop
####################################################################################



###############################
#### calculate tube radius ####
###############################

def nearby_point(n_iteration):
    '''
    This function calculates location of grid points at less than 1 unit above or below current sheet
    Case 1: ## At j, the closest mesh point on the left of sheet is i_left_at_j
            ## At j+1, the closest mesh point on the left of sheet is at least 1 unit greater than i_left_at_j.
    Case 2: ## At j, the closest mesh point on the left of sheet is i_left_at_j
            ## At j+1, the closest mesh point on the left of sheet is still the same i_left_at_j
    '''
    
    if n_iteration == 0:
        print('n_iteration given is not a posive integer')
        return 'Error'
    
    Lower_Right = []
    Upper_Left  = []
    tube_grid_point = np.array([tube_radius(j, n_iteration) for j in range(iv.N_side)])
    
    for j in range(iv.N_side-1):
        i_right = int( tube_grid_point[j+1] )
        i_left  = int( tube_grid_point[j] )
        n_line  = i_right - i_left
        
        if n_line >= 1: 
            ## This is case-1:
            for idx_line in range(1, n_line+1):
                Lower_Right += [ [i_left+idx_line, j] ]     # the mesh grid point below the current sheet 
                Upper_Left  += [ [i_left+idx_line, j+1] ]   # the mesh grid point above the current sheet 
    
    np.save(iv.B_LT_output_dir+'nearby_point/Lower_Right%i.npy'%(n_iteration), Lower_Right)
    np.save(iv.B_LT_output_dir+'nearby_point/Upper_Left%i.npy'%(n_iteration), Upper_Left)
    return None            

################
def J_star(i_ctr, j, n_iteration):
    '''
    Surface-weighted current at the current sheet.
    i_ctr = r_tube/hr
    '''
    if n_iteration <= 0:
        print('n_iteration should be a positive integer')
        return 'Error'
    
    if i_ctr < int(tube_radius(j, n_iteration)) or i_ctr >= int(tube_radius(j, n_iteration))+1:
        return 0

    z_global_bottom_top = np.array([iv.z_bottom, iv.z_bottom + j*iv.hz], dtype=np.float64)
    exp_sh = PM.exponential_factor_scale_height(z_global_bottom_top)
    
    PB_ext_int_diff = B_LT.B_internal_lattice(int(iv.R_star/iv.hr), 0, n_iteration)**2 / (8*np.pi)  # this equals (P_ext-P_int) at the base of tube
    numerator = exp_sh * PB_ext_int_diff                                                 # this equals (P_ext-P_int) at grid height j
    denominator = B_LT.B_internal_lattice(int(i_ctr), j, n_iteration) + B_LT.B_external_lattice(int(i_ctr)+1, j, n_iteration)  # B_int + B_ext at grid height j
    return 2.0 * numerator / denominator

################
def calculate_current_density(n_iteration):
    '''
    This function converts sheet current into volume current.
    '''
    
    if n_iteration <= 0:
        print('n_iteration should be a positive integer')
        return 'Error'
    
    ## execute nearby_point() for the n_th iteration ##
    nearby_point(n_iteration)
    
    tube_grid_point = np.array([tube_radius(j, n_iteration) for j in range(iv.N_side)])
    np.save(iv.B_LT_output_dir+'tube_grid_point/tube_grid_point_%i.npy'%(n_iteration), (tube_grid_point) )
    
    tube_grid_point_left = tube_grid_point.astype(int)
    tube_grid_point_right = tube_grid_point.astype(int)+1
    
    ## Apply forward difference method for sigma_arr; works decently, so no need to use central difference method ##
    sigma_arr = np.arctan( (tube_grid_point[1:] - tube_grid_point[:-1]) * iv.hr/iv.hz )    # unit is [rad]
    sigma_arr = np.append(sigma_arr, 0)  
    
    J_star_arr = np.array([J_star(tube_grid_point[j], j, n_iteration) for j in range(iv.N_side)])    
    l1_grid = tube_grid_point - tube_grid_point_left
    l2_grid = tube_grid_point_right - tube_grid_point
    l1 = l1_grid * iv.hr
    l2 = l2_grid * iv.hr
    J1 = l2/(l1+l2)**2 / np.cos(sigma_arr) * J_star_arr
    J2 = l1/(l1+l2)**2 / np.cos(sigma_arr) * J_star_arr
    return J1, J2, tube_grid_point_left, tube_grid_point_right

################
def calculate_source_cs(n_iteration):
    '''
    This function calculates the source contributions from current sheet
    '''
    if n_iteration == 0:
        return np.zeros(( iv.Nx_matrix*iv.Ny_matrix, 1))
    
    if n_iteration >= 1:
        source = np.zeros(( iv.Nx_matrix*iv.Ny_matrix, 1))
        J1, J2, tube_grid_point_left, tube_grid_point_right = calculate_current_density(n_iteration)
        for j in range(1, iv.N_side):
            i_left = tube_grid_point_left[j]
            i_right = tube_grid_point_right[j]
            i_left = np.where(i_left == 0, 1, i_left)
            i_right = np.where(i_right == iv.N_side, iv.N_side-1, i_right)
            idx_source_left = (i_left-1)*iv.Ny_matrix + j - 1
            idx_source_right = (i_right-1)*iv.Ny_matrix + j - 1
            source[idx_source_left]  = -4*np.pi*(i_left*iv.hr)*J1[j]*iv.hr*iv.hz
            source[idx_source_right] = -4*np.pi*(i_right*iv.hr)*J2[j]*iv.hr*iv.hz
        return source
####################################################################