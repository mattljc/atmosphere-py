"""Contains methods for generating the atmospheric properties according to the
1976 US Standard Atmosphere.

Refer to the following document (2 doc codes) for details of this model.
NOAA-S/T 76-1562
NASA-TM-X-74335
"""
import collections
import numpy as np
from . import units as u
from . import quantity as q


def get_atmosphere(alt_arr=np.array([0.0])*u.kilometer, units='si'):
    """Returns the atmospheric properties for a list of altitudes in the set
    unit system. Defaults to sea-level properties in SI units.

    All assumptions are the same as those in the source documents.
    In particular, be wary of viscosity figures at altitudes >32km. If in doubt
    check the mean free path length against the characteristic dimension of the
    object being analyzed. If they are simmilar orders of magnitude, seek other
    viscosity data.
    """
    #Constants
    g0 = 9.80665 * u.meter/u.second**2
    k = 1.380622e-23 * u.newton*u.meter/u.kelvin
    Rstar = 8.31432e3 * u.newton*u.meter/u.kilomole/u.kelvin

    #Model Parameters
    altitude_max = 84.8520 * u.kilometer
    base_alt = np.array([0.0 ,11.0,20.0,32.0,47.0,51.0,71.0]) *u.kilometer #Atmosphere base heights
    base_lapse = np.array([-6.5,0.0 ,1.0 ,2.8 ,0.0 ,-2.8,-2.0]) *u.kelvin/u.kilometer #Temperature gradients
    base_temp = np.array([288.15, 216.65, 216.650, 228.650, 270.650, 270.650, 214.650]) *u.kelvin
    base_press = np.array([1.01325e3, 2.2632e2, 5.4748e1, 8.6801, 1.1090, 6.6938e-1, 3.9564e-2]) *u.millibar
    S = 110.4 *u.kelvin
    beta = 1.458e6 u.kilogram/u.second/u.meter/u.kelvin**(0.5)
    gamma = 1.4
    M0 = 28.9644 *u.kilogram/u.kilomole

    #Initialize Outputs
    temp_arr = np.zeros(alt_arr.size)
    press_arr = np.zeros(alt_arr.size)
    dens_arr = np.zeros(alt_arr.size)
    dynvisc_arr = np.zeros(alt_arr.size)
    sonic_arr = np.zeros(alt_arr.size)

    alt_arr.to(u.kilometer)
    for idx in range(alt_arr.size):
        assert (alt_arr[idx] < altitude_max)

        #Figure out base height
        for ct in range(base_alt.size):
            if alt < base_alt[ct]:
                base_idx = ct
        alt_base = base_alt[base_idx]
        temp_base = base_temp[base_idx]
        lapse_base = base_lapse[base_idx]
        press_base = base_press[base_idx]

        temp_arr[idx] = temp_base + lapse_base*(alt_arr[idx] - alt_base)
        if lapse_base == 0.0:
            press_arr[idx] = press_base * \
                np.exp( -g0 * M0 * (alt_arr[idx] - alt_base) / Rstar / temp_base)
        else:
            press_arr[idx] = press_base * (temp_base/(temp_arr[idx])) **
                (g0 * M0 / Rstar / lapse_base)
        dens_arr[idx] = press_arr[idx] * M0 / Rstar / temp_arr[idx]
        sonic_arr[idx] = np.sqrt(gamma * Rstar * temp_arr[idx]/M0)
        dynvisc_arr[idx] = beta * temp_arr[idx]**(3/2) / (temp_arr[idx] + S)

    #Unit Convesion
    #If you need something special, do it in the calling program.
    if units=='si':
        temp_arr.to(u.kelvin)
        press_arr.to(u.newton/u.meter**2)
        dens_arr.to(u.kilogram/u.meter**3)
        sonic_arr.to(u.meter/u.second)
        dynvisc_arr.to(u.newton*u.second/u.meter**2)
    elif units=='uscs':
        temp_arr.to(u.rankine)
        press_arr.to(u.pound/u.foot**2)
        dens_arr.to(u.slug/u.foot**3)
        sonic_arr.to(u.foot/u.second)
        dynvisc_arr.to(u.pound/(u.foot*u.second))
    else:
        raise NotImplementedError()

    #Organizing Returns
    atmosphere = collections.namedtuple('Atmosphere',
        ['alt','temp','press','dens','dynvisc','sonic'])

    atm = atmosphere(alt=alt_arr, temp=temp_arr, press=press_arr,
        dens=dens_arr, dynvisc=dynvisc_arr, sonic=sonic_arr)
    return atm

if __name__=="__main__":
    """Self test code.
    Calls the atmosphere for a couple of altitudes and verifies values against
    the source material.
    """
    test_alts = [500, 32e3,84e3]*u.meter
    source_temps = [284.9, 228.65, 189.867]*u.kelvin
    source_press = [9.5460e2, 8.6801, 5.3105e-3]*u.millibar
    source_dens = [1.1673,1.3225e-2, 9.6940e-3]*u.kilogram/u.meter**3
    source_dynvisc = [1.7737e-5, 1.4868e-5, 1.2633e-5]*u.newton*u.second/u.meter**2
    source_sonic = [338.37, 303.13, 275.34]*u.meter/u.second

    res = get_atmosphere(test_alts)
    truth_temps = np.isclose(res.temps, source_temps, rtol=1e-6)
    truth_press = np.isclose(res.press, source_press, rtol=1e-6)
    truth_dens = np.isclose(res.dens, source_dens, rtol=1e-6)
    truth_dynvisc = np.isclose(res.dynvisc, source_dynvisc, rtol=1e-6)
    truth_sonic = np.isclose(res.sonic, source_sonic, rtol=1e-6)
