import numpy as np
from useful import UsefulFunctions

# Obliquity of the ecliptic at J2000.0
EPSILON = np.radians(23.43927944)  # 23° 26' 21.406''

useful = UsefulFunctions()

#                    --- Useful functions ---

    # Getting the Spherical coords
def input_spherical(a):

    print("\n\t--- Enter the coordinates ---")

    r = float(input("Distance (AU):\t"))

    # --- HELIOCENTRIC ---
    # Ecliptic Longitude: 0° to 360°
    # Ecliptic Latitude : -90° to +90°

    if a == 1:

        print("\nEcliptical Longitude (DMS)")

        d_long = float(input("Degrees (0° to 360°): \t"))
        dm_long = float(input("Arcminutes: \t"))
        ds_long = float(input("Arcseconds: \t"))

        useful.verify_dms(360, d_long, dm_long, ds_long)

        print("\nEcliptical Latitude (DMS)")

        d_lat = float(input("Degrees (-90° to +90°): \t"))
        dm_lat = float(input("Arcminutes: \t"))
        ds_lat = float(input("Arcseconds: \t"))

        useful.verify_dms(90, d_lat, dm_lat, ds_lat)

        long_deg = useful.dms_to_degrees(d_long, dm_long, ds_long)
        lat_deg = useful.dms_to_degrees(d_lat, dm_lat, ds_lat)

        long_rad = np.radians(long_deg)
        lat_rad = np.radians(lat_deg)

        return np.array([r, long_rad, lat_rad])

    # --- GEOCENTRIC ---
    # Right Ascension : 0h to 24h
    # Declination     : -90° to +90°

    elif a == 2:

        print("\nRight Ascension (HH MM SS):")

        h = float(input("Hours: \t"))
        m = float(input("Minutes: \t"))
        s = float(input("Seconds: \t"))

        useful.verify_hms(h, m, s)

        print("\nDeclination (DD MM SS)")

        d = float(input("Degrees (use negative if needed): \t"))
        dm = float(input("Arcminutes: \t"))
        ds = float(input("Arcseconds: \t"))

        useful.verify_dms(90, d, dm, ds)

        ra_deg = useful.hms_to_degrees(h, m, s)
        dec_deg = useful.dms_to_degrees(d, dm, ds)

        ra_rad = np.radians(ra_deg)
        dec_rad = np.radians(dec_deg)

        return np.array([r, ra_rad, dec_rad])


#                    --- Coordinate transformations ---

class Coords:

    """
    Coordinate systems:

    Origin:
        1 : Heliocentric
        2 : Geocentric
    Plane:
        1 : Ecliptic
        2 : Equatorial
    Type:
        1 : Cartesian
        2 : Spherical
    """

    def __init__(self, coords, earth_coords,
                 decision_earth,
                 from_origin,
                 from_plane,
                 from_type):

        self.coords = coords
        self.earth_coords = earth_coords
        self.decision_earth = decision_earth
        self.from_origin = from_origin
        self.from_plane = from_plane
        self.from_type = from_type

    #                --- Geometric transformations ---

    def cartesian_to_spherical(self, vector):

        x, y, z = vector
        r = np.linalg.norm(vector)

        lon = np.arctan2(y, x)
        lon = lon % (2*np.pi)

        lat = np.arctan2(z, np.sqrt(x**2 + y**2))

        return np.array([r, lon, lat])

    def spherical_to_cartesian(self, vector):

        r, lon, lat = vector

        x = r * np.cos(lat) * np.cos(lon)
        y = r * np.cos(lat) * np.sin(lon)
        z = r * np.sin(lat)

        return np.array([x, y, z])

    #                --- Rotations ---

    def rotate_x(self, vector, angle):

        c = np.cos(angle)
        s = np.sin(angle)

        R = np.array([
            [1, 0, 0],
            [0, c, -s],
            [0, s, c]
        ])

        return R @ vector

    def ecliptic_to_equatorial(self, vector):
        return self.rotate_x(vector, EPSILON)

    def equatorial_to_ecliptic(self, vector):
        return self.rotate_x(vector, -EPSILON)

    #                --- Origins ---

    def helio_to_geo(self, vector):
        return vector - self.earth_coords

    def geo_to_helio(self, vector):
        return vector + self.earth_coords

    #                --- General change ---

    def general_change(self):

        # Convert to the base coordinates:
        # Cartesian / Ecliptic / Heliocentric

        coords = self.coords.copy()

        if self.from_type == 2:
            coords = self.spherical_to_cartesian(coords)
        if self.from_plane == 2:
            coords = self.equatorial_to_ecliptic(coords)

        base_origin = self.from_origin

        if self.from_origin == 2 and self.decision_earth == 1:
            coords = self.geo_to_helio(coords)
            base_origin = 1

        base = coords.copy()

        results = {}

        origins = {
            1: "Heliocentric",
            2: "Geocentric"
        }
        planes = {
            1: "Ecliptic",
            2: "Equatorial"
        }
        types = {
            1: "Cartesian",
            2: "Spherical"
        }

        # Convert to all coordinate systems

        for o in [1, 2]:
            if self.decision_earth == 2 and o != self.from_origin:
                continue
            for p in [1, 2]:
                for t in [1, 2]:
                    temp = base.copy()
                    # Origin
                    if o != base_origin and self.decision_earth == 1:
                        temp = self.helio_to_geo(temp)
                    # Plane
                    if p == 2:
                        temp = self.ecliptic_to_equatorial(temp)
                    # Type
                    if t == 2:
                        temp = self.cartesian_to_spherical(temp)

                    key = f"{origins[o]} - {planes[p]} - {types[t]}"

                    results[key] = temp

        return results

def manual_input():
    # receive the coordinates from the excecution insted of the terminal
    print("\t--- CURRENT COORDINATES ---")
    print("\nWhat coordinates do you have?")

    # ORIGIN
    a = int(input("\nHeliocentric (1)    Geocentric (2): \t"))
    useful.verify_1or2(a)

    # PLANE
    b = int(input("\nEcliptic (1)    Equatorial (2): \t"))
    useful.verify_1or2(b)

    # TYPE
    c = int(input("\nCartesian (1)    Spherical (2): \t"))
    useful.verify_1or2(c)

    # Current coordinates

        # Object coords
    if c == 2:

        initial_coords = input_spherical(a)
        useful.verify_coords(initial_coords)

    else:
        initial_coords = np.array(list(map(float,input(
            "\nEnter the 3 OBJECT coordinates separated by space:\n").split())))
        useful.verify_coords(initial_coords)

        # Earth coords

    earth_coords = np.array([0.0, 0.0, 0.0])

    e = int(input(
        "\nDo you have the Earth coordinates?"
        "\nYES (1)    NO (2)\n"))
    useful.verify_1or2(e)

    if e == 1:
        earth_coords = np.array(list(map(float,input(
                    "\nEnter the 3 EARTH coordinates separated by space:\n").split())))
        useful.verify_coords(earth_coords)

    return initial_coords, earth_coords, a, b, c, e


#                    --- Main program ---

def main():
    
    try:

        initial_coords, earth_coords, a, b, c, e = manual_input()

        # Coordinate transformation

        converter = Coords(
            coords=initial_coords,
            earth_coords=earth_coords,
            decision_earth=e,
            from_origin=a,
            from_plane=b,
            from_type=c
        )

        results = converter.general_change()

        # Print results

        print("\n\t--- ALL COORDINATE SYSTEMS ---\n")

        for k, v in results.items():

            if "Spherical" in k and "Geocentric" in k:

                r, ra, dec = v
                h, hm, hs = useful.radians_to_hms(ra)
                sign, d, dm, ds = useful.radians_to_dms(dec)
                sign_str = "-" if sign < 0 else "+"
                print(k)
                print(f"  Distance: {r:.8f} AU")
                print(f"  RA : {h}h {hm}m {hs:.3f}s")
                print(f"  Dec: {sign_str}{d}° {dm}' {ds:.3f}''\n")

            elif "Spherical" in k and "Heliocentric" in k:

                r, lon, lat = v
                sign_lon, d_lon, dm_lon, ds_lon = useful.radians_to_dms(lon)
                sign_lat, d_lat, dm_lat, ds_lat = useful.radians_to_dms(lat)
                sign_lon_str = "-" if sign_lon < 0 else "+"
                sign_lat_str = "-" if sign_lat < 0 else "+"
                print(k)
                print(f"  Distance: {r:.8f} AU")
                print(
                    f"  Longitude: "
                    f"{sign_lon_str}{d_lon}° {dm_lon}' {ds_lon:.3f}''"
                )
                print(
                    f"  Latitude : "
                    f"{sign_lat_str}{d_lat}° {dm_lat}' {ds_lat:.3f}''\n"
                )

            else:
                print(k)
                print(v, "\n")

    except Exception as error:
        print(f"\nERROR: {error}")


if __name__ == "__main__":
    main()