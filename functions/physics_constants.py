'''
Sun, in CGS units
'''
R_sun = 6.957e10                 # Radius of the Sun, [cm]
M_sun = 1.989e33                 # Mass of the Sun, [gram]
AU    = 1.496e13                 # 1 AU, [cm]
sq_Rsun_AU = (R_sun/AU)**2       # square of R_sun/AU

'''
Constants, in CGS units
'''
c_light = 29979245800.         # speed of light, [cm/s]
m_p = 1.6726219e-24            # proton mass, [gram]
m_H = 1.6735328e-24            # hydrogen mass, [gram]
m_e = 9.1093837e-28            # electron mass, [gram]
e_charge = 4.80320427e-10      # electric charge, [StatC]
kb = 1.380649e-16              # Boltzmann constant, [erg/K]
G_N = 6.6743e-8                # Gravitational constant, [dyn cm^2 / g^2]
m_u = 1.660539e-24             # Atomic mass unit (AMU) [gram]
two_pi_c = 188365156730.8853   # (2) times (numpy.pi) times (speed of light), [cm/s]
A_H = 1.0                      # mass number of hydrogen
A_He = 4.0                     # mass number of helium

'''
Unit conversion
'''
Tesla2Gauss = 1e4              ## Convert Tesla to Gauss; one Tesla is 10,000 Gauss
Gauss2Tesla = 1e-4             ## Convert Gauss to Tesla; one Gauss is 1e-4 Tesla
GeV2erg = 1.6021773e-3         ## Convert GeV to erg; one GeV is 1.602e-3 erg
erg2GeV = 6.2415065e2
cm2km = 1e-5                   ## Convert cm to km; one cm is 1e-5 km
mbarn2cmsq = 1e-27             ## convert micro Barn to cm^2; 1 Barn = 1e-28 m^2


'''
pp cross-section
delta-functional method
'''
K_pi = 0.17
n_tilde = 1.23409
m_p_GeV =  0.9382721           ## proton rest mass * c^2, unit in [GeV]
m_pi_GeV = 0.1349768           ## pion rest mass*c**2, unit is GeV
E_th_GeV = m_p_GeV + 2.0*m_pi_GeV + m_pi_GeV**2/(2.0*m_p_GeV)     ## E_th_GeV = 1.2179 GeV
Fr_delta_reduced = 2.0 * n_tilde / K_pi