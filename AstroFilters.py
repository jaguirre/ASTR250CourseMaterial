import numpy as np
#import matplotlib.pyplot as plt
import h5py
from astropy import units as u #, constants as c
from scipy.interpolate import CubicSpline

# This is going to be a problem ...
path = '../.'

filter_names = ['U','B', 'V', 'R', 'I']

def ReadJohnsonCousins():
    jc_filter_data = {}
    f = h5py.File("johnson_cousins.hdf5", "r")
    for filter_name in filter_names:
        jc_filter_data[filter_name] = {}
        jc_filter_data[filter_name]['wvl'] = f[filter_name]['wvl'][:].copy()
        jc_filter_data[filter_name]['filt'] = f[filter_name]['filt'][:].copy()
    f.close()
    return jc_filter_data

jc_filter_data = ReadJohnsonCousins()

def InterpolateFilter(filtername):
    wvl = jc_filter_data[filtername]['wvl'].copy()
    filt = jc_filter_data[filtername]['filt'].copy()

    wvl /= 10. # nm
    # Need to pad with zeros outside the defined range
    wvl_min = 200.
    wvl_max = 1100.
    nwvl = 1000
    iwvl = np.linspace(wvl_min,wvl_max,num=nwvl)
    ifilt = np.zeros(nwvl)
    wh = np.logical_and(iwvl>wvl.min(),iwvl<wvl.max())
    ifilt[wh] = CubicSpline(wvl,filt)(iwvl[wh]) #np.interp(iwvl[wh],wvl,filt)
    synthfilt = CubicSpline(iwvl,ifilt)
    
    return synthfilt

# Build the filter dictionary (soon to object)
filters = {}
for filtername in filter_names:
    filters[filtername] = {}
    filters[filtername]['shape'] = InterpolateFilter(filtername)

wave = np.linspace(200,1100,num=1000)

for fn in filter_names:
    dwvl = np.trapz(filters[fn]['shape'](wave),x=wave)*u.nm
    wvl_mean = np.trapz(wave*filters[fn]['shape'](wave),x=wave)/dwvl*u.nm**2
    filters[fn]['wvl_mean'] = wvl_mean
    filters[fn]['dwvl'] = dwvl

zero_unit = u.erg/u.s/u.cm**2/u.Angstrom
filters['U']['zero'] = 4.19e-9 * zero_unit
filters['B']['zero'] = 6.60e-9 * zero_unit
filters['V']['zero'] = 3.55e-9 * zero_unit
filters['R']['zero'] = 2.25e-9 * zero_unit
filters['I']['zero'] = 1.22e-9 * zero_unit
