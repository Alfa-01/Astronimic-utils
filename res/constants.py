"""Module containing constants for calculations"""

import numpy as np
import numpy.typing as npt


# Constants
MU: float = 398600.4418
EARTH_RADIUS: float = 6371.0
OMEGA_EARTH: float = 7.2921159e-5

# Earth's ellipsoid parameters (WGS84)
a: float = 6378137.0
f: float = 1 / 298.257223563
b: float = a * (1 - f)
e: float = 1 - (b**2 / a**2)

# Orbital elements
SEMIMAJOR_AXIS: float = 42_164
ECCENTRICITY: float = 0.24
INCLINATION: float = 1.11
RAAN: float = 2.27
PERIGEE_ARGUMENT: float = 4.71

# Args
PRECISION_FOR_KEPLER_EQUATION: float = 1e-6
AMOUNT_OF_POINTS: int = 100

# numpy types
NDArrayFloat = npt.NDArray[np.float64]
