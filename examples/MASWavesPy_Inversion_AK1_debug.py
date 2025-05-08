import maswavespy as sw
from maswavespy import inversion

import pandas as pd
import numpy as np


# Path to sample data
# Experimental dispersion curve (stored as a .txt file in the folder Data)
filename_dc = '/Users/tourei@mines.edu/coding/sits/cryoseismic_imaging/maswavespy/data/AK_dc_some_picks.txt' 
# Initial values for soil model parameters (stored as a .csv file in the folder Data)
filename_initial = '/Users/tourei@mines.edu/coding/sits/cryoseismic_imaging/maswavespy/data/AK_initial.csv' 


# Import experimental dispersion curves
wavelengths = []
c_mean = []; c_low = []; c_up = []
with open(filename_dc, 'r') as file_dc:
    next(file_dc) # Skip the header
    for value in file_dc.readlines():
        wavelengths.append(float(value.split()[0]))
        c_mean.append(float(value.split()[1]))
        c_low.append(float(value.split()[2]))
        c_up.append(float(value.split()[3]))
wavelengths = np.array(wavelengths, dtype='float64')
c_mean = np.array(c_mean, dtype='float64'); 
c_low = np.array(c_low, dtype='float64')
c_up = np.array(c_up, dtype='float64')
    
# Import initial soil model parameters
initial_parameters = pd.read_csv(filename_initial)
h = np.array(initial_parameters['h [m]'].values[0:-1], dtype='float64')
n = int(len(h))
Vs = np.array(initial_parameters['Vs [m/s]'].values, dtype='float64')
rho = np.array(initial_parameters['rho [kg/m3]'].values, dtype='float64')
Vp = []
n_unsat = 0; nu = None
for item in range(len(initial_parameters['saturated/unsaturated'].values)):
    if initial_parameters['saturated/unsaturated'].values[item] == 'unsat':
        nu = initial_parameters['nu [-]'].values[item]
        Vp.append(np.sqrt((2*(1-nu))/(1-2*nu))*Vs[item])
        n_unsat = n_unsat + 1
    else:
        Vp.append(initial_parameters['Vp [m/s]'].values[item])
Vp = np.array(Vp, dtype='float64')

# Print message to user
print('The sample dispersion curve has been imported.')
print('The initial soil model parameters have been imported.')


# Initialize an inversion object.    
site = 'AK'
profile = 'P1'
inv_TestSite = inversion.InvertDC(site, profile, c_mean, c_low, c_up, wavelengths)

# Print message to user
print('An inversion (InvertDC) object has been initialized.') 



# Initialize the inversion routine. The inversion is conducted using 
# Monte Carlo sampling as described in Olafsdottir et al. (2020).

# Range for testing phase velocity
c_min = 200; c_max = 2500; c_step = 0.1; delta_c = 3
c_test = {'min' : c_min, 
          'max' : c_max,
          'step' : c_step,
          'delta_c' : delta_c}

# Initial model parameters
initial = {'n' : n,
           'n_unsat' : n_unsat,
           'alpha' : Vp, # alpha unsat
           'nu_unsat' : 0.18,
           'alpha_sat' : 4000, 
           'beta' : Vs,
           'rho' : rho,
           'h' : h,
           'reversals' : 0}

# Inversion algorithm settings. See further in Olafsdottir et al. (2020).
settings = {'run' : 20,
            'bs' : 5,
            'bh' : 10,
            'N_max' : 1000}        

# View the initial shear wave velocity profile.
# Compute the associated dispersion curve and show relative to the experimental
# data. The misfit value is printed to the screen.
max_depth = 140
inv_TestSite.view_initial(initial, max_depth, c_test, col='crimson', DC_yaxis='linear', 
                 fig=None, ax=None, figwidth=16, figheight=12, return_ct=False)

# Print message to user
print('The initial estimate of the Vs profile and the corresponding theoretical DC have been plotted.')




# Start the inversion analysis (optimization) process.
print('Inversion initiated.')
inv_TestSite.mc_inversion(c_test, initial, settings)

