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
    --reference_time: Reference time (t_r) in Julian Date (JD) (e.g., 2461174.5 for May 14,2026 0h UTC)
"""

from coordinates import Coords
from useful import UsefulFunctions, K, RELATIVE_MASSES, ORBITAL_ELEMENTS
import numpy as np 
import argparse

# Reference time "May 14,2026 0h UTC (JD 2461174.5)"
REFERENCE_TIME = 2461174.5

#                --- Solution of the problem ---

class TwoBodyProblem:

    """
    This class contains methods to solve the two-body problem, which describes the motion of two
    celestial bodies under their mutual gravitational attraction. It includes methods to calculate the
    position and velocity vectors of the object, as well as its orbital elements.
    """
    
    def __init__(self, object_name, verbose=True, orbital_elements=None, relative_mass=0, reference_time=0):
        
        self.verbose = verbose
        self.object_name = object_name
        self.a = orbital_elements["a"]
        self.e = orbital_elements["e"]
        self.i = np.radians(orbital_elements["i"])
        self.Omega = np.radians(orbital_elements["Omega"])
        self.omega = np.radians(orbital_elements["omega"])
        self.M_r = np.radians(orbital_elements["M_r"])
        self.relative_mass = relative_mass
        self.t_r = reference_time

        mu = (K**2) * (1 + self.relative_mass)  # Gravitational parameter
        self.n = np.sqrt(mu/ (self.a**3))

    def mean_anomaly(self, t):
        # Calculate the mean anomaly (M) at time t using the formula: M = M_r + n * (t - reference_time)
            # n is the mean motion, calculated as n = K * (180/pi) * sqrt((1 + relative_mass) / a^3)
        
        M_rad = self.M_r + self.n * (t - self.t_r)
        M_deg = np.degrees(M_rad) % 360  # Convert to degrees and ensure it's between 0 and 360

        if self.verbose:
            print(f"Mean Anomaly: (M)= {M_deg:.8f}° \n")
            
        return M_deg

    def eliptic_eccentric_anomaly_newton(self, M_deg, tolerance=1e-8):
        # Solves the Kepler equation for the elipctic case using the Newton-Raphson method 
        # to find the eccentric anomaly (E) given the mean anomaly (M) and eccentricity (e).
        M_rad = np.radians(M_deg)
        E_rad = M_rad  # Initial value
        
        while True:
            # f(E) = E - e*sin(E) - M
            f_E = E_rad - self.e * np.sin(E_rad) - M_rad
            # Derivada: f'(E) = 1 - e*cos(E)
            f_prime_E = 1 - self.e * np.cos(E_rad)
            
            # Newton-Raphson
            E_new = E_rad - (f_E / f_prime_E)
            if abs(E_new - E_rad) < tolerance:
                break
            E_rad = E_new
            
        E_deg = np.degrees(E_rad)
        if self.verbose:
            print(f"Eccentric Anomaly: (E) = {E_deg:.8f}° \n")
        return E_deg

    def position_vector(self, E):
        # Calculate the position vector (r) using the formula: r = a * (1 - e * cos(E))
        r = self.a * (1 - self.e * np.cos(np.radians(E)))
        if self.verbose:
            print(f"Position Vector: (r)= {r:.8f} au \n")
        return r

    def true_anomaly(self, E):
        # Calculate the true anomaly (theta) using the formula: 
        # theta = 2 * atan2(sqrt((1 + e)/(1-e)) * tan2(E/2))
        theta = 2 * np.arctan2(np.sqrt((1 + self.e)/(1-self.e)) * np.tan(np.radians(E) / 2), 1)
        if self.verbose:
            print(f"True Anomaly: (theta)= {np.degrees(theta):.8f}° \n")
        return theta

    def position_vector_cartesian(self, r, theta):
        # Calculate the position vector in Cartesian coordinates (x, y, z) using the formulas:
        # x = r * (cos(Omega) * cos(theta + omega) - sin(Omega) * sin(theta + omega) * cos(i))
        # y = r * (sin(Omega) * cos(theta + omega) + cos(Omega) * sin(theta + omega) * cos(i))
        # z = r * (sin(theta + omega) * sin(i))
        x = r * (np.cos(self.Omega) * np.cos(theta + self.omega) - 
                np.sin(self.Omega) * np.sin(theta + self.omega) * np.cos(self.i))
        
        y = r * (np.sin(self.Omega) * np.cos(theta + self.omega) + 
                np.cos(self.Omega) * np.sin(theta + self.omega) * np.cos(self.i))
        
        z = r * (np.sin(theta + self.omega) * np.sin(self.i))
        
        if self.verbose:
            print(f"Position Vector in Heliocentric-Ecliptical-Cartesian Coordinates: (->r)=[{x:.8f}, {y:.8f}, {z:.8f}] au \n")
        return x, y, z

    def velocity_vector_cartesian(self, r, theta):
        # Calculate the velocity vector in Cartesian coordinates (vx, vy, vz) using the formulas:
        # vx = S_11 (r_dot * cos(theta) - r * theta_dot * sin(theta)) + S_12 (r_dot * sin(theta) + r * theta_dot * cos(theta)
        # vy = S_21 (r_dot * cos(theta) - r * theta_dot * sin(theta)) + S_22 (r_dot * sin(theta) + r * theta_dot * cos(theta))
        # vz = S_31 (r_dot * cos(theta) - r * theta_dot * sin(theta)) + S_32 (r_dot * sin(theta) + r * theta_dot * cos(theta))
        # where S_ij are the elements of the transformation matrix from perifocal to ecliptic coordinates, r_dot is the radial velocity, and theta_dot is the angular velocity
        # r_dot = (K * sqrt(1 + relative_mass) * e* sin(theta))/(sqrt(a * (1 - e^2)))
        # theta_dot = ((K * sqrt(1 + relative_mass) * (sqrt(a * (1 - e^2)) / (r^2)))
        r_dot = (K * np.sqrt(1 + self.relative_mass) * self.e * np.sin(theta)) / (np.sqrt(self.a * (1 - self.e**2)))
        theta_dot = (K * np.sqrt(1 + self.relative_mass) * (np.sqrt(self.a * (1 - self.e**2)) / (r**2)))

        S_11 = np.cos(self.Omega) * np.cos(self.omega) - np.sin(self.Omega) * np.sin(self.omega) * np.cos(self.i)
        S_12 = -np.cos(self.Omega) * np.sin(self.omega) - np.sin(self.Omega) * np.cos(self.omega) * np.cos(self.i)
        S_21 = np.sin(self.Omega) * np.cos(self.omega) + np.cos(self.Omega) * np.sin(self.omega) * np.cos(self.i)
        S_22 = -np.sin(self.Omega) * np.sin(self.omega) + np.cos(self.Omega) * np.cos(self.omega) * np.cos(self.i)
        S_31 = np.sin(self.omega) * np.sin(self.i)
        S_32 = np.cos(self.omega) * np.sin(self.i)

        vx = S_11 * (r_dot * np.cos(theta) - r * theta_dot * np.sin(theta)) + S_12 * (r_dot * np.sin(theta) + r * theta_dot * np.cos(theta))
        vy = S_21 * (r_dot * np.cos(theta) - r * theta_dot * np.sin(theta)) + S_22 * (r_dot * np.sin(theta) + r * theta_dot * np.cos(theta))
        vz = S_31 * (r_dot * np.cos(theta) - r * theta_dot * np.sin(theta)) + S_32 * (r_dot * np.sin(theta) + r * theta_dot * np.cos(theta))

        if self.verbose:
            print(f"Velocity Vector in Heliocentric-Ecliptical-Cartesian Coordinates: (->v)= [{vx:.8f}, {vy:.8f}, {vz:.8f}] au/dia \n")
            print(f"Velocity Magnitude: (v)= {r_dot:.8f} au/dia \n")
        return vx, vy, vz, r_dot

    def general_solution (self, t):
        M = self.mean_anomaly(t)
        E = self.eliptic_eccentric_anomaly_newton(M)
        r = self.position_vector(E)
        theta = self.true_anomaly(E)
        x, y, z = self.position_vector_cartesian(r, theta)
        vx, vy, vz, v = self.velocity_vector_cartesian(r, theta)

        return M, E, r, theta, x, y, z, vx, vy, vz, v

#                --- Useful functions (Manual Input) ---

def manual_input():

    print("\t--- MANUAL INPUT ---")
    object_name = input("\nEnter the name of the object in the solar system (e.g., Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune):\n")

    if object_name not in ORBITAL_ELEMENTS:
        print(f"Object '{object_name}' is not in the database.")
        a = float(input("\nEnter the semi-major axis (a) (in astronomical units, au):\n"))
        e = float(input("\nEnter the eccentricity (e) (dimensionless):\n"))
        i = float(input("\nEnter the inclination (i) (in degrees):\n"))
        Omega = float(input("\nEnter the longitude of the ascending node (Omega) (in degrees):\n"))
        omega = float(input("\nEnter the argument of periapsis (omega) (in degrees):\n"))
        M_r = float(input("\nEnter the mean anomaly of reference (M_r) (in degrees):\n"))
        reference_time = float(input("\nEnter the reference time (t_r) in Julian Date (JD) (e.g., 2461174.5 for May 14,2026 0h UTC):\n"))

    else:
        a = ORBITAL_ELEMENTS[object_name]["a"]
        e = ORBITAL_ELEMENTS[object_name]["e"]
        i = ORBITAL_ELEMENTS[object_name]["i"]
        Omega = ORBITAL_ELEMENTS[object_name]["Omega"]
        omega = ORBITAL_ELEMENTS[object_name]["omega"]
        M_r = ORBITAL_ELEMENTS[object_name]["M_r"]

    return object_name, a, e, i, Omega, omega, M_r, reference_time

#              --- Useful functions (converters) ---

useful = UsefulFunctions()

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
    parser.add_argument('--reference_time', type=float, help='Reference time (t_r) in Julian Date (JD) (e.g., 2461174.5 for May 14,2026 0h UTC)')
    args = parser.parse_args()

    print("\t########## SOLUTION OF THE TWO-BODY PROBLEM ########## \n")

    #   --- Initialize the orbital elements ---
    if args.object is None:
        object_name, a, e, i, Omega, omega, M_r , reference_time = manual_input()
    else:
        object_name = args.object
        a = args.a if args.a is not None else ORBITAL_ELEMENTS[object_name]["a"]
        e = args.e if args.e is not None else ORBITAL_ELEMENTS[object_name]["e"]
        i = args.i if args.i is not None else ORBITAL_ELEMENTS[object_name]["i"]
        Omega = args.Omega if args.Omega is not None else ORBITAL_ELEMENTS[object_name]["Omega"]
        omega = args.omega if args.omega is not None else ORBITAL_ELEMENTS[object_name]["omega"]
        M_r = args.M_r if args.M_r is not None else ORBITAL_ELEMENTS[object_name]["M_r"]
        reference_time = args.reference_time if args.reference_time is not None else REFERENCE_TIME 
    
    if args.time is None:
        time = float(input("\nEnter the time in Julian Date (e.g., 2461174.5 for May 14,2026 0h UTC):\n"))
    else:        
        time = args.time

    #  --- Solution of the problem --- 

    print(f"--- Object: {object_name}, Time: {time} (JD) ---\n")

    # Solution of the problem for the object at the given time
    object_elements = {"a": a, "e": e, "i": i, "Omega": Omega, "omega": omega, "M_r": M_r}

    object_solution = TwoBodyProblem(
                    object_name=object_name, 
                    orbital_elements=object_elements, 
                    relative_mass=RELATIVE_MASSES.get(object_name, 0), 
                    reference_time=reference_time)

    M, E, r, theta, x, y, z, vx, vy, vz, v = object_solution.general_solution(time)

    print(f"--- Object: Earth, Time: {time} (JD) ---\n")

    # Solution of the problem for the Earth at the given time 
    # (to calculate the position of the object with respect to the Earth)
    earth_solution = TwoBodyProblem(
                    object_name="Earth", 
                    orbital_elements=ORBITAL_ELEMENTS["Earth"], 
                    relative_mass=RELATIVE_MASSES["Earth+Moon"], 
                    reference_time=reference_time)
    
    M_earth, E_earth, r_earth, theta_earth, x_earth, y_earth, z_earth, vx_earth, vy_earth, vz_earth , v_earth = earth_solution.general_solution(time)

    print(f"--- Position of {object_name} with respect to the Earth ---\n")

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
    h, hm, hs = useful.radians_to_hms(r_geo_equ_sph[1])
    sign, d, dm, ds = useful.radians_to_dms(r_geo_equ_sph[2])
    string_str = "-" if sign < 0 else "+"
    print(f"Position Vector in Geocentric-Equatorial-Spherical Coordinates: (->r)=[{r:.8f} AU, {h}h {hm}m {hs:.3f}s, {string_str}{d}°{dm}' {ds:.3f}''] \n")


if __name__ == "__main__":
    main()