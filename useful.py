import numpy as np

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