import os

import numpy as np
import h5py
import matplotlib.pyplot as plt

def coolingrate(tables, density, He_mass_frac=0.258, metallicity=0.25):
    elements = [
        x for x in tables.keys() if x not in [
            'Header', 'Metal_free', 'Solar', 'Total_Metals']
    ]
    
    metal_free_He_mass_frac = list(
        tables['Metal_free']['Helium_mass_fraction_bins']
    )
    metal_free_ne_nh = list(
        tables['Metal_free']['Electron_density_over_n_h']
    )
    
    He_mass_frac_index = metal_free_He_mass_frac.index(np.float32(He_mass_frac))
    abundance_index = list(
        tables[elements[0]]['Hydrogen_density_bins']).index(np.float32(density)
    )

    lambda_metal_free = tables['Metal_free']\
                              ['Net_Cooling']\
                              [He_mass_frac_index, :, abundance_index]
    
    T = tables['Metal_free']['Temperature_bins']
    
    sum_metals = np.zeros(352)
    for e in elements:
        lambda_i_solar = tables[e]['Net_Cooling'][:,abundance_index]
        
        metal_free_ne_nh = tables['Metal_free']\
                                 ['Electron_density_over_n_h']\
                                 [He_mass_frac_index, :,abundance_index]
        solar_ne_nh = tables['Solar']\
                            ['Electron_density_over_n_h']\
                            [:, abundance_index]
        
        sum_metals += lambda_i_solar*(metal_free_ne_nh/solar_ne_nh)*metallicity
        
    coolingrate = lambda_metal_free + sum_metals
    
    return T, coolingrate

def interpolate(y1, y2, x1, x2, x, log=False):
    if log:
        y1 = np.log10(y1)
        y2 = np.log10(y2)
        x1 = np.log10(x1)
        x2 = np.log10(x2)
        x = np.log10(x)
    y = (y2 - y1) / (x2 - x1) * (x - x1) + y1
    
    if log:
        y = 10**y
    return y

if __name__ == '__main__':
    datadir = 'data/CoolingTables/'
    print('1a')
    f1 = h5py.File("../exercise1/data/CoolingTables/z_3.017.hdf5", 'r')
    f2 = h5py.File("../exercise1/data/CoolingTables/z_2.829.hdf5", 'r')
    
    for i in [1, 1e-2, 1e-4, 1e-6]:
        T, L_z2829 = coolingrate(
                         tables=f1, 
                         density=i, 
                         He_mass_frac=0.258, 
                         metallicity=0.25
                     )
        T, L_z3017 = coolingrate(
                         tables=f2, 
                         density=i, 
                         He_mass_frac=0.258, 
                         metallicity=0.25
                     )
        
        z2829 = 2.829
        z3017 = 3.017
        z = 3
        
        L_z3000 = interpolate(L_z2829, L_z3017, z2829, z3017, z)
        
        plt.loglog(T, L_z3000, label=r'$\rho={}$'.format(i))

    plt.ylabel(
        r'$\Lambda \  [\mathrm{erg} \  \mathrm{cm}^{-3} \  \mathrm{s}^{-1}]$'
    )
    plt.xlabel(r'$T \ [\mathrm{K}]$')
    plt.axis(xmin=5e3, xmax=1e9)
    plt.legend()
    plt.tight_layout()
    plt.savefig('plots/1a_coolingrates.png')

    f1.close()
    f2.close()