import numpy as np
from xarray import apply_ufunc

# Basic functions
def calc_esat_from_tref(tref):
    '''Calculate the saturation vapor pressure [mbar] from the absolute 
    temperature [K].
    
    Parameters
    ----------
    tref : np.array
        Absolute temperate [K]
        
    Returns
    ----------
    esat : np.array
        Saturation vapor pressure [mbar]
    
    '''
    return np.exp(-2991.2729/tref**2 - 6017.0128/tref 
              + 18.87643854 - 0.028354721 * tref
              + 1.7838301E-5*tref**2 - 8.4150417E-10*tref**3
              + 4.4412542E-13*tref**4 + 2.858487 * np.log(tref))/100

def calc_wsat_from_esat_p(esat,p):
    '''Calculate the saturation mixing ratio [g/kg] from the saturation 
    vapor pressure [mbar] and the pressure [mbar].'''
    return 621.97*esat/(p-esat)

def calc_rh_from_sh_wsat(sh,wsat):
    '''Calculate the relative humidity [%] from the specific humidity [kg/kg]
    and the saturation mixing ratio [g/kg].'''
    return 100*sh/((1-sh)*wsat/1000)

def calc_w_from_sh(sh):
    '''Calculate mixing ratio [g/kg] from specific humidity [kg/kg].
    Return in g/kg.
    Ref: https://www.e-education.psu.edu/meteo300/node/519'''
    return 1000*sh/(1-sh)

def calc_w_from_wsat_rh(wsat,rh):
    '''Calculate the mixing ratio [g/kg] from the saturation mixing ratio [g/kg]
    and the relative humidity [%].'''
    return rh/100*wsat

def calc_tl_from_tref_rh(tref,rh):
    '''Calculate the lifting condensation temperature [K] from the absolute 
    temperature [K] and the relative humidity [%].'''
    return 1/(1/(tref-55)-np.log(rh/100)/2840)+55

def calc_theta_from_w_tref_tl_p(w,tref,tl,p):
    '''Calculate the potential temperature equivalent [K] of the lifting 
    condensation temperature from the absolute temperature [K], the lifting 
    condensation temperature [K], and the pressure [mbar].
    
    This is the temperature a parcel of air would reach if it were continued 
    to be adiabatically lifted to condense all water, and then lowered dry 
    adiabatically to 1,000 mbar.'''
    return (tref*(1000/p)**(0.2854*(1-0.28E-3*w))
            *np.exp((3.376/tl-0.00254)*w*(1+0.81E-3*w)))

def calc_wbt_from_theta(theta):
    '''Calculate the wet bulb temperature [C] from the potential temperature 
    equivalent of the lifting condensation temperature [K].'''
    return 45.114-51.489*(theta/273.15)**-3.504

def calc_wbt_from_tref_rh(tref, rh):
    '''Calculate wet bulb temperature from absolute temperature and relative humidity,
    according to the formula of Stull (2011).
    Ref: https://journals.ametsoc.org/doi/full/10.1175/JAMC-D-11-0143.1'''
    tref = tref-273.15
    term1_stull = tref * np.arctan(0.151977 * np.sqrt(rh + 8.313659))
    term2_stull = np.arctan(tref + rh)
    term3_stull = np.arctan(rh - 1.676331)
    term4_stull = 0.00391838 * np.power(rh, 1.5) * np.arctan(0.023101 * rh)
    term5_stull = 4.686035
    return term1_stull + term2_stull - term3_stull + term4_stull - term5_stull

# Higher-level functions
def _calc_wbt_from_tref_sh_p(tref,sh,p,method='Davies-Jones'):
    '''Calculate wet bulb temperature from absolute temperature, specific 
    humidity, and pressure.'''
    if method=='Davies-Jones':
        esat = calc_esat_from_tref(tref)
        wsat = calc_wsat_from_esat_p(esat,p)
        rh = calc_rh_from_sh_wsat(sh,wsat)
        w = calc_w_from_wsat_rh(wsat,rh)
        tl = calc_tl_from_tref_rh(tref,rh)
        theta = calc_theta_from_w_tref_tl_p(w,tref,tl,p)
        wbt = calc_wbt_from_theta(theta)
    elif method=='Stull':
        esat = calc_esat_from_tref(tref)
        wsat = calc_wsat_from_esat_p(esat,p)
        rh = calc_rh_from_sh_wsat(sh,wsat)
        wbt = calc_wbt_from_tref_rh(tref, rh)
        return wbt
    
def _calc_rh_from_tref_sh_p(tref,sh,p):
    '''Calculate relative humidity from absolute temperature, specific 
    humidity, and pressure.'''
    esat = calc_esat_from_tref(tref)
    wsat = calc_wsat_from_esat_p(esat,p)
    rh = calc_rh_from_sh_wsat(sh,wsat)
    return rh

def _calc_wbt_from_tref_rh_p(tref,rh,p):
    '''Calculate wet bulb temperature from absolute temperature [K], relative 
    humidity [%], and pressure [mbar].'''
    esat = calc_esat_from_tref(tref)
    wsat = calc_wsat_from_esat_p(esat,p)
    w = calc_w_from_wsat_rh(wsat,rh)
    tl = calc_tl_from_tref_rh(tref,rh)
    theta = calc_theta_from_w_tref_tl_p(w,tref,tl,p)
    wbt = calc_wbt_from_theta(theta)
    return wbt

# xarray wrappers
def calc_wbt_from_tref_sh_p(temperature_at_2m,specific_humidity,pressure,method):
    '''Calculate wet bulb temperature from absolute temperature, specific 
    humidity, and pressure.'''
    return apply_ufunc(_calc_wbt_from_tref_sh_p,
                         temperature_at_2m,specific_humidity,pressure,method,
                         dask='parallelized',output_dtypes=[temperature_at_2m.dtype])

def calc_rh_from_tref_sh_p(temperature_at_2m, specific_humidity,pressure):
    '''Calculate relative humidity from absolute temperature, specific 
    humidity, and pressure.'''
    return apply_ufunc(_calc_rh_from_tref_sh_p,
                         temperature_at_2m,specific_humidity,pressure,
                         dask='parallelized',output_dtypes=[temperature_at_2m.dtype])  