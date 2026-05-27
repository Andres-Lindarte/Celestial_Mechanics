import numpy as np

class UsefulFunctions:

    """This class contains useful functions for celestial mechanics calculations.
    It includes methods for verifying input values, 
    converting between different coordinate systems, 
    and performing various calculations related to orbital mechanics.
    
    """

    def __init__(self,
                vector=None,        # to verify_coords
                h=None, m=None, s=None,      # to verify_hms, hms_to_degrees
                n=None,            # to verify_1or2, verify_dms
                d=None, #m=None, s=None,     # to verify_dms, dms_to_degrees
                rad=None,          # to radians_to_dms, radians_to_hms
                ):
        
        self.vector = vector
        self.h = h
        self.m = m
        self.s = s
        self.n = n
        self.d = d
        self.rad = rad  

            # Verifications
    def verify_1or2(self):
        if self.n not in (1, 2):
            raise ValueError("Error: Number must be 1 or 2")

    def verify_coords(self):
        if self.vector.size != 3:
            raise ValueError("Error: Coordinates must be a size 3 vector.")

    def verify_dms(self):  # n is the maximum value, i.e., 90 or 360

        if self.n == 90:
            if abs(self.d) > 90:
                raise ValueError("Degrees must be between -90 and 90.")
        elif self.n == 360:
            if not (0 <= self.d < 360):
                raise ValueError("Degrees must be between 0 and 360.")
        if not (0 <= self.m < 60):
            raise ValueError("Arcminutes must be between 0 and 60.")
        if not (0 <= self.s < 60):
            raise ValueError("Arcseconds must be between 0 and 60.")

        return True


    def verify_hms(self):

        if not (0 <= self.h < 24):
            raise ValueError("Hours must be between 0 and 24.")
        if not (0 <= self.m < 60):
            raise ValueError("Minutes must be between 0 and 60.")
        if not (0 <= self.s < 60):
            raise ValueError("Seconds must be between 0 and 60.")

        return True


        # Conversions

            # From Degrees-Minutes-Seconds (DMS) to Degrees
    def dms_to_degrees(self):

        sign = -1 if self.d < 0 else 1
        d = abs(self.d)

        return sign * (d + self.m/60 + self.s/3600)


            # From Hours-Minutes-Seconds (HMS) to Degrees
    def hms_to_degrees(self):

        return 15 * (self.h + self.m/60 + self.s/3600)


            # From Radians to DMS
    def radians_to_dms(self):

        deg = np.degrees(self.rad)
        sign = -1 if deg < 0 else 1
        deg = abs(deg)

        d = int(deg)
        m = int((deg - d) * 60)
        s = (deg - d - m/60) * 3600

        return sign, d, m, s


            # From Radians to HMS
    def radians_to_hms(self):

        rad = self.rad % (2*np.pi)

        total_hours = np.degrees(rad) / 15

        h = int(total_hours)
        m = int((total_hours - h) * 60)
        s = (total_hours - h - m/60) * 3600

        return h, m, s