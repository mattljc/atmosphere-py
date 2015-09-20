"""Contains methods for generating the atmospheric properties according to the
1976 US Standard Atmosphere.

Refer to the following document (2 doc codes) for details of this model.
NOAA-S/T 76-1562
NASA-TM-X-74335
"""
import numpy as np
from . import units as u
from . import quantity as q


"""Returns the atmospheric properties for a list of altitudes in the set unit
system.

Defaults to sea-level properties in SI units. The notation follows that used in
the original source document.
"""
def getAtmosphere(H=np.array([0.0]), units='si'):
    #Constants
    g0 = 9.80665 * u.meter/u.second**2
    k = 1.380622e-23 * u.newton*u.meter/u.kelvin
    Rstar = 8.31432e3 * u.newton*u.meter/u.kilomole/u.kelvin

    #Model Parameters
    altitude_max = 84.8520 * u.kilometer
    Hb = np.array([11.0,20.0,32.0,47.0,51.0,71.0]) *u.kilometer #Atmosphere base heights
    Lb = np.array([-6.5,0.0,1.0,2.8,0.0,-2.8,-0.0]) *u.kelvin/u.kilometer #Temperature gradients
    S = 110 *u.kelvin
    beta = 1.458e6 u.kilogram/u.second/u.meter/u.kelvin**(0.5)
    gamma = 1.4
    P0 = 1.013250e5 *u.newton/u.meter**2
    T0 = 288.15 *u.kelvin
    M0 = 28.9644 *u.kilogram/u.kilomole

    for Hi in H:
        assert (Hi < altitude_max)
        #Temperature
        ct = 0
        T = T0
        while (H<Hb[ct]):
            T = T + L[ct]*(H-H[ct])
