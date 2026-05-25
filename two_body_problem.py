"""
Usage:
    python two_body_problem.py --object "Mars" --time 2461174.5

Description:
    This program calculates the position of an object (especially a planet) in the solar system 
    at a given time (in Julian Date) using the two-body problem. It takes into account the orbital 
    elements of the object and the Earth to compute the position in various coordinate systems 
    (heliocentric-ecliptical-cartesian, geocentric-ecliptical-cartesian, geocentric-equatorial-cartesian, 
    and geocentric-equatorial-spherical coordinates). 

    # Manuel input:
    You can give the orbital elements manually or use the predefined values for the planets in the solar system.
    
    # Command-line arguments (CLI):
    You can use the command-line arguments:
    --object: Name of the object (e.g., Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune)
    --time: Time in Julian Date (e.g., 2461174.5 for May 14,2026 0h UTC)
    --a: Semi-major axis (in astronomical units, au)
    --e: Eccentricity (dimensionless)
    --i: Inclination (in degrees)
    --Omega: Longitude of the ascending node (in degrees)
    --omega: Argument of periapsis (in degrees)
    --M_r: Mean anomaly of reference (in degrees)
"""

from coordenates import Coords
import numpy as np 
import argparse

# Gaussian gravitational constant in au^(3/2)/dia
K = 0.01720209908  
# Inclination of the ecliptic from the celestial equator
EPSILON = np.radians(23.43927944)  # 23° 26' 21.406'' 
# Reference time "May 14,2026 0h UTC (JD 2461174.5)"
reference_time = 2461174.5

#                  --- Useful directories ---     

# Relative masses of the planets with respect to the Sun
relative_masses = {
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
orbital_elements = {
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

#                --- Useful functions (Solution of the problem) ---

def mean_excentric_anomaly(M_r, t, a, object_name, relative_mass=0):
    # Calculate the mean anomaly (M) at time t using the formula: M = M_r + n * (t - reference_time)
        # n is the mean motion, calculated as n = K * (180/pi) * sqrt((1 + relative_mass) / a^3)
    if object_name in relative_masses:
        relative_mass = relative_masses[object_name]
    n = K * (180/np.pi) * np.sqrt((1 + relative_mass) / a**3)
    
    M = M_r + n * (t - reference_time)
    print(f"Mean Anomaly of {object_name} in {t} (JD): (M)= {M:.8f}° \n")
    return M

def eliptic_eccentric_anomaly_newton(M_deg, e, tolerance=1e-8):
    # Solves the Kepler equation for the elipctic case using the Newton-Raphson method 
    # to find the eccentric anomaly (E) given the mean anomaly (M) and eccentricity (e).
    M_rad = np.radians(M_deg)
    E_rad = M_rad  # Initial value
    
    while True:
        # f(E) = E - e*sin(E) - M
        f_E = E_rad - e * np.sin(E_rad) - M_rad
        # Derivada: f'(E) = 1 - e*cos(E)
        f_prime_E = 1 - e * np.cos(E_rad)
        
        # Newton-Raphson
        E_new = E_rad - (f_E / f_prime_E)
        if abs(E_new - E_rad) < tolerance:
            break
        E_rad = E_new
        
    E_deg = np.degrees(E_rad)
    print(f"Eccentric Anomaly: (E) = {E_deg:.8f}° \n")
    return E_deg


def position_vector(a, e, E):
    # Calculate the position vector (r) using the formula: r = a * (1 - e * cos(E))
    r = a * (1 - e * np.cos(np.radians(E)))
    print(f"Position Vector: (r)= {r:.8f} au \n")
    return r

def true_anomaly(E, e):
    # Calculate the true anomaly (theta) using the formula: 
    # theta = 2 * atan2(sqrt((1 + e)/(1-e)) * tan2(E/2))
    theta = 2 * np.arctan(np.sqrt((1 + e)/(1-e)) * np.tan(np.radians(E) / 2))
    print(f"True Anomaly: (theta)= {np.degrees(theta):.8f}° \n")
    return theta

def position_vector_cartesian(r, theta, i, Omega, omega):
    # Calculate the position vector in Cartesian coordinates (x, y, z) using the formulas:
    # x = r * (cos(Omega) * cos(theta + omega) - sin(Omega) * sin(theta + omega) * cos(i))
    # y = r * (sin(Omega) * cos(theta + omega) + cos(Omega) * sin(theta + omega) * cos(i))
    # z = r * (sin(theta + omega) * sin(i))
    x = r * (np.cos(np.radians(Omega)) * np.cos(theta + np.radians(omega)) - 
             np.sin(np.radians(Omega)) * np.sin(theta + np.radians(omega)) * np.cos(np.radians(i)))
    
    y = r * (np.sin(np.radians(Omega)) * np.cos(theta + np.radians(omega)) + 
             np.cos(np.radians(Omega)) * np.sin(theta + np.radians(omega)) * np.cos(np.radians(i)))
    
    z = r * (np.sin(theta + np.radians(omega)) * np.sin(np.radians(i)))
    
    print(f"Position Vector in Heliocentric-Ecliptical-Cartesian Coordinates: (->r)=[{x:.8f}, {y:.8f}, {z:.8f}] au \n")
    return x, y, z

#                --- Useful functions (Manual Input) ---

def manual_input():

    print("\t--- MANUAL INPUT ---")
    object_name = input("\nEnter the name of the object in the solar system (e.g., Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune):\n")

    if object_name not in orbital_elements:
        print(f"Object '{object_name}' is not in the database.")
        a = float(input("\nEnter the semi-major axis (a) (in astronomical units, au):\n"))
        e = float(input("\nEnter the eccentricity (e) (dimensionless):\n"))
        i = float(input("\nEnter the inclination (i) (in degrees):\n"))
        Omega = float(input("\nEnter the longitude of the ascending node (Omega) (in degrees):\n"))
        omega = float(input("\nEnter the argument of periapsis (omega) (in degrees):\n"))
        M_r = float(input("\nEnter the mean anomaly of reference (M_r) (in degrees):\n"))

    else:
        a = orbital_elements[object_name]["a"]
        e = orbital_elements[object_name]["e"]
        i = orbital_elements[object_name]["i"]
        Omega = orbital_elements[object_name]["Omega"]
        omega = orbital_elements[object_name]["omega"]
        M_r = orbital_elements[object_name]["M_r"]

    return object_name, a, e, i, Omega, omega, M_r

#              --- Useful functions (converters) ---

        # From Radians to DMS
def radians_to_dms(rad):

    deg = np.degrees(rad)
    sign = -1 if deg < 0 else 1
    deg = abs(deg)

    d = int(deg)
    m = int((deg - d) * 60)
    s = (deg - d - m/60) * 3600

    return sign, d, m, s


        # From Radians to HMS
def radians_to_hms(rad):

    rad = rad % (2*np.pi)

    total_hours = np.degrees(rad) / 15

    h = int(total_hours)
    m = int((total_hours - h) * 60)
    s = (total_hours - h - m/60) * 3600

    return h, m, s

#             --- Main function ---

def main ():
    #   --- Command-line arguments (CLI) OPTIONAL ---
    parser = argparse.ArgumentParser(description='Calculate the position of an object (especially a planet) in the solar system.')
    parser.add_argument('--object', type=str, help='Name of the object (e.g., Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune)')
    parser.add_argument('--time', type=float, help='Time in Julian Date (e.g., 2461174.5 for May 14,2026 0h UTC)')
    parser.add_argument('--a', type=float, help='Semi-major axis (in astronomical units, au)')
    parser.add_argument('--e', type=float, help='Eccentricity (dimensionless)')
    parser.add_argument('--i', type=float, help='Inclination (in degrees)')
    parser.add_argument('--Omega', type=float, help='Longitude of the ascending node (in degrees)')
    parser.add_argument('--omega', type=float, help='Argument of periapsis (in degrees)')
    parser.add_argument('--M_r', type=float, help='Mean anomaly of reference (in degrees)')
    args = parser.parse_args()

    #   --- Initialize the orbital elements ---
    if args.object is None:
        object_name, a, e, i, Omega, omega, M_r = manual_input()
    else:
        object_name = args.object
        a = args.a if args.a is not None else orbital_elements[object_name]["a"]
        e = args.e if args.e is not None else orbital_elements[object_name]["e"]
        i = args.i if args.i is not None else orbital_elements[object_name]["i"]
        Omega = args.Omega if args.Omega is not None else orbital_elements[object_name]["Omega"]
        omega = args.omega if args.omega is not None else orbital_elements[object_name]["omega"]
        M_r = args.M_r if args.M_r is not None else orbital_elements[object_name]["M_r"]

    if args.time is None:
        time = float(input("\nEnter the time in Julian Date (e.g., 2461174.5 for May 14,2026 0h UTC):\n"))
    else:        
        time = args.time

    #  --- Solution of the problem --- 
    
    # Solution of the problem for the object at the given time
    M = mean_excentric_anomaly(M_r, time, a, object_name)
    E = eliptic_eccentric_anomaly_newton(M, e)
    r = position_vector(a, e, E)
    theta = true_anomaly(E, e)
    x, y, z = position_vector_cartesian(r, theta, i, Omega, omega)

    # Solution of the problem for the Earth at the given time 
    # (to calculate the position of the object with respect to the Earth)
    M_earth = mean_excentric_anomaly(orbital_elements["Earth"]["M_r"], time, orbital_elements["Earth"]["a"], "Earth")
    E_earth = eliptic_eccentric_anomaly_newton(M_earth, orbital_elements["Earth"]["e"])
    r_earth = position_vector(orbital_elements["Earth"]["a"], orbital_elements["Earth"]["e"], E_earth)
    theta_earth = true_anomaly(E_earth, orbital_elements["Earth"]["e"])
    x_earth, y_earth, z_earth = position_vector_cartesian(r_earth, theta_earth, orbital_elements["Earth"]["i"], orbital_elements["Earth"]["Omega"], orbital_elements["Earth"]["omega"])

    # Position of the object with respect to the Earth
    converter = Coords(
        coords=0,
        earth_coords=np.array([x_earth, y_earth, z_earth]),
        decision_earth=0,
        from_origin=0,
        from_plane=0,
        from_type=0
    )

    r_geo = converter.helio_to_geo(np.array([x, y, z]))
    print(f"Position Vector in Geocentric-Ecliptical-Cartesian Coordinates: (->r)=[{r_geo[0]:.8f}, {r_geo[1]:.8f}, {r_geo[2]:.8f}] au \n")

    r_geo_equ = converter.ecliptic_to_equatorial(r_geo)
    print(f"Position Vector in Geocentric-Equatorial-Cartesian Coordinates: (->r)=[{r_geo_equ[0]:.8f}, {r_geo_equ[1]:.8f}, {r_geo_equ[2]:.8f}] au \n")

    r_geo_equ_sph = converter.cartesian_to_spherical(r_geo_equ)
    r = r_geo_equ_sph[0]
    h, hm, hs = radians_to_hms(r_geo_equ_sph[1])
    sign, d, dm, ds = radians_to_dms(r_geo_equ_sph[2])
    string_str = "-" if sign < 0 else "+"
    print(f"Position Vector in Geocentric-Equatorial-Spherical Coordinates: (->r)=[{r:.8f} AU, {h}h {hm}m {hs:.3f}s, {string_str}{d}°{dm}' {ds:.3f}''] \n")


if __name__ == "__main__":
    main()