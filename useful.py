"""
Description:
    This code is not excecutable by its onw.
    This code gathers useful constants:
    - Gaussian gravitational constant (K)
    - Relative masses of the planets with respect to the Sun
    - Orbital elements of the planets for May 14,2026 0h UTC

    It also provides a class with useful functions for celestial mechanics calculations, 
    including methods for:
    - Verifying input values: 
        - verify_1or2: Ensures the input is either 1 or 2.
        - verify_coords: Ensures the input is a size 3 vector.
        - verify_dms: Ensures the input degrees, arcminutes, and arcseconds are within valid ranges.
        - verify_hms: Ensures the input hours, minutes, and seconds are within valid ranges.
    - Converting between different coordinate systems:
        - dms_to_degrees: Converts degrees, arcminutes, and arcseconds to decimal degrees.
        - hms_to_degrees: Converts hours, minutes, and seconds to decimal degrees.
        - radians_to_dms: Converts radians to degrees, arcminutes, and arcseconds.
"""

import numpy as np

#                   --- Useful quantities ---

# Obliquity of the ecliptic at J2000.0
EPSILON = np.radians(23.43927944)  # 23° 26' 21.406''

# Gaussian gravitational constant in au^(3/2)/dia
K = 0.01720209908  

# Relative masses of the planets with respect to the Sun
RELATIVE_MASSES = {
    "Mercury": 1/6023600,
    "Venus": 1/408523.71,
    "Earth+Moon": 1/328900.56,
    "Mars": 1/3098708.0,
    "Jupiter": 1/1047.3486,
    "Saturn": 1/3497.898,
    "Uranus": 1/22902.98,
    "Neptune": 1/19412.24,
}

# Orbital elements of the planets for May 14,2026 0h UTC (JD 2461174.5)
"""
a : Semi-major axis (in astronomical units, au)
e : Eccentricity (dimensionless)
i : Inclination (in degrees)
Omega : Longitude of the ascending node (in degrees)
omega : Argument of periapsis (in degrees)
M_r : Mean anomaly of reference (in degrees)
"""
ORBITAL_ELEMENTS = {
    "Mercury": {"a": 0.38709817, 
                "e": 0.20563319,
                "i": 7.00338683,
                "Omega": 48.29789723,
                "omega": 29.20193758,
                "M_r": 341.92845078
                },
    "Venus": {"a": 0.72333150, 
              "e": 0.00679081,
              "i": 3.39437486,
              "Omega": 76.60678210,
              "omega": 54.87676238,
              "M_r": 358.21130004
              },
    # Earth + Moon 
    "Earth": {"a": 1.00074332,
              "e": 0.01639114,
              "i": 0.00284898,
              "Omega": 146.48865232,
              "omega": 318.95788506,
              "M_r": 125.91233118
                },
    "Mars": {"a": 1.52362081,
             "e": 0.09343917,
             "i": 1.84749141,
             "Omega": 49.48112365,
             "omega": 286.61454328,
             "M_r": 25.52866375
                },
    "Jupiter": {"a": 5.20182186,
                "e": 0.04825346,
                "i": 1.30379948,
                "Omega": 100.52125276,
                "omega": 273.38083997,
                "M_r": 100.59906410
                },
    "Saturn": {"a": 9.54167409,
               "e": 0.05515038,
               "i": 2.48823401,
               "Omega": 113.65630093,
               "omega": 338.54518281,
               "M_r": 280.13785163
                },
    "Uranus": {"a": 19.27448311,
               "e": 0.04721564,
               "i": 0.77202058,
               "Omega": 74.02932312,
               "omega": 91.83592007,
               "M_r": 255.18851944
                },
    "Neptune": {"a": 30.09568563,
                "e": 0.01050175,
                "i": 1.77110623,
                "Omega": 131.82116887,
                "omega": 277.76383938,
                "M_r": 312.22914269
                },
}

#                   --- Useful functions ---

class UsefulFunctions:

    """This class contains useful functions for celestial mechanics calculations.
    It includes methods for verifying input values, 
    converting between different coordinate systems, 
    and performing various calculations related to orbital mechanics.
    
    """

    @staticmethod
    def verify_1or2(n):
        if n not in (1, 2):
            raise ValueError("Error: Number must be 1 or 2")

    
    @staticmethod
    def verify_coords(vector):
        if vector.size != 3:
            raise ValueError("Error: Coordinates must be a size 3 vector.")

    @staticmethod
    def verify_dms(n, d, m, s):  # n is the maximum value, i.e., 90 or 360

        if n == 90:
            if abs(d) > 90:
                raise ValueError("Degrees must be between -90 and 90.")
        elif n == 360:
            if not (0 <= d < 360):
                raise ValueError("Degrees must be between 0 and 360.")
        if not (0 <= m < 60):
            raise ValueError("Arcminutes must be between 0 and 60.")
        if not (0 <= s < 60):
            raise ValueError("Arcseconds must be between 0 and 60.")

        return True


    @staticmethod
    def verify_hms(h, m, s):

        if not (0 <= h < 24):
            raise ValueError("Hours must be between 0 and 24.")
        if not (0 <= m < 60):
            raise ValueError("Minutes must be between 0 and 60.")
        if not (0 <= s < 60):
            raise ValueError("Seconds must be between 0 and 60.")

        return True


        # Conversions

            # From Degrees-Minutes-Seconds (DMS) to Degrees
    @staticmethod
    def dms_to_degrees(d, m, s):

        sign = -1 if d < 0 else 1
        d = abs(d)

        return sign * (d + m/60 + s/3600)


            # From Hours-Minutes-Seconds (HMS) to Degrees
    @staticmethod
    def hms_to_degrees(h, m, s):

        return 15 * (h + m/60 + s/3600)


            # From Radians to DMS
    @staticmethod
    def radians_to_dms(rad):

        deg = np.degrees(rad)
        sign = -1 if deg < 0 else 1
        deg = abs(deg)

        d = int(deg)
        m = int((deg - d) * 60)
        s = (deg - d - m/60) * 3600

        return sign, d, m, s


            # From Radians to HMS
    @staticmethod
    def radians_to_hms(rad):

        rad = rad % (2*np.pi)

        total_hours = np.degrees(rad) / 15

        h = int(total_hours)
        m = int((total_hours - h) * 60)
        s = (total_hours - h - m/60) * 3600

        return h, m, s