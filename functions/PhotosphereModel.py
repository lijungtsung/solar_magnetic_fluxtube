import numpy as np
from functions import physics_constants as const
from scipy import integrate

#### load data ####
data_atmospheric = np.loadtxt('data_input/data_SolarAtmosphere/solar_atmosphere_profile.txt')

#### reverse array; make them start from the smallest quantity ####
z_SL = data_atmospheric[:,0][::-1]    # height, from small to large
T_SL = data_atmospheric[:,1][::-1]    # temperature, from small to large
P_gas_SL = data_atmospheric[:,2][::-1]    # gas pressure, from small to large
rho_tot_SL = data_atmospheric[:,3][::-1]   # total mass density, from small to large
nH_SL = data_atmospheric[:,4][::-1]   # hydrogen number density, from small to large
MU_SL = data_atmospheric[:,5][::-1]   # mean molecular weight MU, from small to large

if np.all(np.diff(z_SL) >= 0) != True:
    print('Error: z_SL is not monotonically increasing array')

#### max height in the data ####
data_z_max = max(data_atmospheric[:,0])
data_z_min = min(data_atmospheric[:,0])

#### rewrite nH and rho functions, by using np.interp, which is faster ####
def n_as_z(z):
    ## hydrogen number density, as function of z
    '''
    Numpy interpolation is faster than scipy interp1d;
    Note that numpy interp only accept monotonically increasing sample points.
    '''
    return np.interp(z, z_SL, nH_SL)

def T_as_z(z):
    ## temperature, as function of z
    '''
    Numpy interpolation is faster than scipy interp1d;
    Note that numpy interp only accept monotonically increasing sample points.
    '''
    return np.interp(z, z_SL, T_SL)

def rho_tot_as_z(z):
    ## total mass density, as function of z
    '''
    Numpy interpolation is faster than scipy interp1d;
    Note that numpy interp only accept monotonically increasing sample points.
    '''
    return np.interp(z, z_SL, rho_tot_SL)

def P_gas_as_z(z):
    ## gas pressure, as function of z
    '''
    Numpy interpolation is faster than scipy interp1d;
    Note that numpy interp only accept monotonically increasing sample points.
    '''
    return np.interp(z, z_SL, P_gas_SL)

def MU_as_z(z):
    ## mean molecular weight MU, as function of z
    '''
    Numpy interpolation is faster than scipy interp1d;
    Note that numpy interp only accept monotonically increasing sample points.
    '''
    return np.interp(z, z_SL, MU_SL)


#### define functions ####
def H_scale_height(z):
    '''
    Calculate the scale height at height z
    Return unit is [cm]
    '''
    
    if np.ndim(z)==0:
        z = np.array([z])
    if np.ndim(z)!=0 and np.ndim(z)!=1:
        print('input z is not an array or scalar')
        return 'input z is not an array or scalar'
    
    if max(z) >= data_z_max:
        print('maximum input z is out of data range')
        return 'maximum input z is out of data range'
    if min(z) <= data_z_min:
        print('minimum input z is out of data range')
        return 'minimum input z is out of data range'
    
    return (const.R_sun)**2 * const.kb * T_as_z(z) / ( const.G_N * const.M_sun * const.m_u * MU_as_z(z) )

def exponential_factor_scale_height(z_start_end):
    
    if np.ndim(z_start_end)!=1:
        print('input z_start_end format is not dim=1')
        return 'input z_start_end format is not dim=1'
    
    z_lower, z_higher = z_start_end
    
    if z_lower > z_higher:
        print('z_lower must be smaller than z_higher')
        return 'z_lower must be smaller than z_higher'
    
    elif z_lower == z_higher:
        return 1.0
    
    else:
        z_integrated = np.linspace(z_lower, z_higher, 1001)
        exponent = integrate.simps(1./H_scale_height(z_integrated), z_integrated, axis=-1)
        return np.exp(-exponent)