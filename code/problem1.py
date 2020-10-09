import os

import h5py
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt

def bisect(x_range, x_i):
    """ Find the indices of elements in x_range that enclose x_i. If x_i outside 
        of x_range, the begin or end index is returned.

    Args:
        x_range (list): Sorted ascending list.
        x_i (float, int): Float or int to find encloding values for

    Raises:
        RuntimeError: No edges found.
    """
    if x_i <= x_range[0]:
        return 0
    if x_i >= x_range[-1]:
        return len(x_range)-1
    
    for i in range(len(x_range)-1):
        if x_i > x_range[i] and x_i < x_range[i+1]:
            return (i, i+1)
    raise RuntimeError("No edges found.")

def coolingrate(tables, density, He_mass_frac=0.258, metallicity=0.25):
    """ Calculate the cooling rate from information in cooling tables. 
        
        The density and He_mass_frac should be in the Hydrogen_density_bins and 
        Metal_free.Helium_mass_fraction_bins, since a direct match is looked up.

    Args:
        tables ([type]): [description]
        density ([type]): [description]
        He_mass_frac (float, optional): [description]. Defaults to 0.258.
        metallicity (float, optional): [description]. Defaults to 0.25.

    Returns:
        List with temperatures and list with cooling rates as function of 
        temperature.
    """

    # list of all elements in table
    elements = [
        x for x in tables.keys() if x not in [
            'Header', 'Metal_free', 'Solar', 'Total_Metals']
    ]
    
    # Helium mass fraction bins of only hydrogen and helium
    metal_free_He_mass_frac = list(
        tables['Metal_free']['Helium_mass_fraction_bins']
    )

    # Ratio of free electron and hydrogen number densities (n_e/n_H) as a 
    # function of /Metal_free/Hydrogen_density_bins, 
    # /Metal_free/Temperature_bins, and helium abundance.
    metal_free_ne_nh = list(
        tables['Metal_free']['Electron_density_over_n_h']
    )
    
    # Find index of He mass frac, ssuming direct match in table.
    He_mass_frac_index = metal_free_He_mass_frac.index(np.float32(He_mass_frac))
    # Find index of density bin, ssuming direct match in table.
    abundance_index = list(
        tables[elements[0]]['Hydrogen_density_bins']).index(np.float32(density)
    )

    # Ratio of free electron and hydrogen number densities (n_e/n_H) as a 
    # function of temperature at given denisty and helium mass fraction.
    metal_free_ne_nh = list(
        tables['Metal_free']\
              ['Electron_density_over_n_h']\
              [He_mass_frac_index, :,abundance_index]
    )

    # Ratio of free electron and hydrogen number densities (n_e/n_H) as a 
    # function of temperature at given density for solar abundance.
    solar_ne_nh = tables['Solar']\
                        ['Electron_density_over_n_h']\
                        [:, abundance_index]

    # Normalized net cooling rate for only hydrogen and helium as a function of 
    # /Metal_free/Hydrogen_density_bins, /Metal_free/Temperature_bins, 
    # and helium abundance. 
    lambda_metal_free = tables['Metal_free']\
                              ['Net_Cooling']\
                              [He_mass_frac_index, :, abundance_index]
    
    # Temperature bins (T [K])
    T = tables['Metal_free']['Temperature_bins']
    
    sum_metals = np.zeros(352)
    for e in elements:
        # Normalized, net cooling rate e as a function of temperature for 
        # Hydrogen density 'density'.
        lambda_i_solar = tables[e]['Net_Cooling'][:,abundance_index]
        
        # Calulating the cooling contribution of element e and addind to total.
        sum_metals += lambda_i_solar*(metal_free_ne_nh/solar_ne_nh)*metallicity
        
    coolingrate = lambda_metal_free + sum_metals
    
    return T, coolingrate

def interpolate(y1, y2, x1, x2, x, log=False):
    """ Linear interpolator. 

    Args:
        y1 (int/float): y value of left interval boundry.
        y2 (int/float): y value of right interval boundry.
        x1 (int/flout): x value of left interval boundry.
        x2 (int/float): x value of right interval boundry.
        x (int/float): x value for which to calculate interpolation.
        log (bool, optional): If true, interpolate linear in log-space. Defaults
             to False.

    Returns:
        float: interpolated y(x)
    """
    if log:
        # transform values to log space
        y1 = np.log10(y1)
        y2 = np.log10(y2)
        x1 = np.log10(x1)
        x2 = np.log10(x2)
        x = np.log10(x)

    # Calculate y(x)
    y = (y2 - y1) / (x2 - x1) * (x - x1) + y1
    
    if log:
        # transform result back
        y = 10**y
    return y

if __name__ == '__main__':
    datadir = 'data/CoolingTables/'
    
    # 1a
    f1 = h5py.File(os.path.join(datadir, "z_3.017.hdf5"), 'r')
    f2 = h5py.File(os.path.join(datadir, "z_2.829.hdf5"), 'r')
    
    # loop over denisties
    for i in [1, 1e-2, 1e-4, 1e-6]:
        # Calculate cooling rate, as function of temperature, at z=2.829.
        T, L_z2829 = coolingrate(
                         tables=f1, 
                         density=i, 
                         He_mass_frac=0.258, 
                         metallicity=0.25
                     )
        # Calculate cooling rate, as function of temperature, at z=3.017.
        T, L_z3017 = coolingrate(
                         tables=f2, 
                         density=i, 
                         He_mass_frac=0.258, 
                         metallicity=0.25
                     )
        
        z2829 = 2.829
        z3017 = 3.017
        z = 3
        
        # interpolate cooling rate for z=3.00.
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
    plt.close()

    f1.close()
    f2.close()

    # 1b
    n = 1e-4 # number density
    metallicity = 0.5 # fraction of solar metallicity
    z_true = sorted([
        float(x[2:7]) for x in os.listdir(datadir) 
            if not any(i in x for i in ['nocompton', 'collis', 'photodis']) 
            and x.endswith('.hdf5')
    ]) # List with true z values in dataset

    # set range for new z
    z_min, z_max, z_step = 0, 9, 100
    z_range = np.linspace(z_min, z_max, z_step)

    L = []
    print("Generating plots for video ...")
    for z_i in tqdm(z_range): 
        
        if z_i in z_true:
            # if z_i in tables, use that insted of interpolating.
            tablename = 'z_{:.3f}.hdf5'.format(z_i)
            tables = h5py.File(os.path.join(datadir, tablename), 'r')
            T, L_i = coolingrate(
                tables=tables, 
                density=n, 
                He_mass_frac=0.258, 
                metallicity=metallicity
            )
            L.append(L_i)
        else:
            # Find neigboring x values
            z_idx = bisect(z_true, z_i)
            # handle edge cases
            if isinstance(z_idx, int) and z_idx == len(z_true)-1:
                z_idx = (z_idx-1, z_idx)

            z_lower, z_upper = z_true[z_idx[0]], z_true[z_idx[1]]
            tablename_lower =  'z_{:.3f}.hdf5'.format(z_lower)
            tablename_upper =  'z_{:.3f}.hdf5'.format(z_upper)
            
            tables_lower = h5py.File(os.path.join(datadir, tablename_lower),'r')
            tables_upper = h5py.File(os.path.join(datadir, tablename_upper),'r')
                    
            
            # Calculate lower cooling rate
            T, L_lower = coolingrate(
                tables=tables_lower, 
                density=n, 
                He_mass_frac=0.258, 
                metallicity=metallicity
            )
            
            # Calculate upper cooling rate
            _, L_upper = coolingrate(
                tables=tables_upper, 
                density=n, 
                He_mass_frac=0.258, 
                metallicity=metallicity
            )

            # interpolate cooling rate for z_i
            L_i = interpolate(L_lower, L_upper, z_lower, z_upper, z_i)
            
            L.append(L_i)
            
        
        plt.loglog(T, L_i, label=r'$z={:.3f}$'.format(z_i))
        plt.text(
            x=1e2, 
            y=1e-20, 
            s=r'$z={:.3f}$'.format(z_i)\
              + '\n' \
              + r'$\rho={:.4f}'.format(n)\
              + r'\ \mathrm{cm}^{-3}$'
        )
        plt.ylabel(
            r'$\Lambda \  [\mathrm{erg} \  \mathrm{cm}^{-3} \  \mathrm{s}^{-1}]$'
        )
        plt.xlabel(r'$T \ [\mathrm{K}]$')
        plt.axis(ymin=1e-27, ymax=1e-19)
        plt.savefig('plots/coolingrate_z{}.png'.format(z_i), dpi=200)
        plt.close()