from find_vectors import FindVectors
import numpy as np
import argparse

# Gaussian gravitational constant in au^(3/2)/dia
K = 0.01720209908  

#                    --- Useful directories ---

relative_masses = {
    # Relative masses of the planets with respect to the Sun
    "Mercury": 1/6023600,
    "Venus": 1/408523.71,
    "Earth+Moon": 1/328900.56,
    "Mars": 1/3098708.0,
    "Jupiter": 1/1047.3486,
    "Saturn": 1/3497.898,
    "Uranus": 1/22902.98,
    "Neptune": 1/19412.24,
}

#                    --- Useful functions ---

find_vectors = FindVectors()

class OrbitalObjects:

    """This class contains methods to calculate important orbital elements of celestial objects,
    such as:
    - Semi-major axis                   (a)
    - Eccentricity                      (e) 
    - Inclination                       (i)
    - Longitude of the ascending node   (Omega)
    - Argument of periapsis             (omega)
    - Mean anomaly                      (M)
    - Proper time of periapsis passage  (t_0)
    
    """
    def position_vector(self, x, y, z):
        return find_vectors.position_vector(x, y, z) 
    # It returns: position, r_magnitude

    def velocity_vector(self, x_speed, y_speed, z_speed):
        return find_vectors.velocity_vector(x_speed, y_speed, z_speed)
    # It returns: velocity, v_magnitude
    
    def angular_momentum_vector(self, position, velocity):
        return find_vectors.angular_momentum_vector(position, velocity)
    # It returns: angular_momentum, h_magnitude
    
    def eccentricity_vector(self, planet, position, velocity, r_magnitude, angular_momentum):
        return find_vectors.eccentricity_vector(planet, position, velocity, r_magnitude, angular_momentum)
    # It returns: e_vector, e_magnitude

    @staticmethod
    def semi_major_axis(planet, v_magnitude, r_magnitude):
        # Calculate the semi-major axis (a) using the formula: a = r / (2-Q), where Q = r * v^2 / mu

        Q = r_magnitude * v_magnitude**2 / K**2 / (1 + relative_masses[planet])  # Q = r * v^2 / mu
        a = r_magnitude / (2 - Q)

        print(f"Semi-major axis: (a)= {a:.8f} ua \n")
        return a

    def rate_of_change_of_position(self, r_magnitude, position, velocity):
        return find_vectors.rate_of_change_of_position(r_magnitude, position, velocity)
    # It returns: dr_dt

    @staticmethod
    def mean_anomaly(a, e_magnitude, r_magnitude, dr_dt):

        E = np.arccos((a - r_magnitude) / (a * e_magnitude))  # Eccentric anomaly
        M = E - e_magnitude * np.sin(E)  # Mean anomaly
        M = np.degrees(M) % 360  # Convert to degrees and ensure it's between 0 and 360

        print(f"Mean anomaly: (M)= {M:.8f} ° \n")
        return M
    
    @staticmethod
    def inclination(h_magnitude, angular_momentum):
        # Calculate the inclination (i) using the formula: i = arccos(h_z / h)
        h_z = angular_momentum[2]
        i = np.arccos(h_z / h_magnitude)
        
        print(f"Inclination: (i)= {np.degrees(i):.8f} ° \n")
        return i
    
    @staticmethod
    def longitude_of_ascending_node(angular_momentum):
        # Calculate the longitude of the ascending node (Omega) using the formula: Omega = arctan2(h_x, -h_y)
        h_x = angular_momentum[0]
        h_y = angular_momentum[1]
        Omega = np.arctan2(h_x, -h_y)
        Omega = np.degrees(Omega) % 360  # Convert to degrees and ensure it's between 0 and 360
        
        print(f"Longitude of the ascending node: (Omega)= {Omega:.8f} ° \n")
        return Omega
    
    @staticmethod
    def argument_of_periapsis(e_vector, angular_momentum, h_magnitude):
        # Calculate the argument of periapsis (omega) using the formula: omega = arctan2(e_z*h_magnitude, -e_x*h_y + e_y*h_x)
        e_x = e_vector[0]
        e_y = e_vector[1]
        e_z = e_vector[2]
        h_x = angular_momentum[0]
        h_y = angular_momentum[1]
        omega = np.arctan2(e_z * h_magnitude, -e_x * h_y + e_y * h_x)
        omega = np.degrees(omega) % 360  # Convert to degrees and ensure it's between 0 and 360
        
        print(f"Argument of periapsis: (omega)= {omega:.8f} ° \n")
        return omega
    
    @staticmethod
    def proper_time_of_periapsis_passage(t_r, M, planet, a):
        # Calculate the proper time of periapsis passage (t_0) using the formula: t_0 = t_r -M/n, where n is the mean motion (n = sqrt(mu/a^3))
        mu = K**2 * (1 + relative_masses[planet])  # Standard gravitational parameter for the planet
        n = np.sqrt(mu / a**3)  # Mean motion
        t_0 = t_r -M / n

        print(f"Proper time of periapsis passage: (t_0)= {t_0:.8f} JD \n")
        return t_0
    
#                --- Manual input ---

def manual_input():

    print("\t--- MANUAL INPUT ---")
    planet = input("Object name (e.g., Mercury, Earth+Moon): \n")

    if planet not in relative_masses:
        print(f"The object '{planet}' is not in the database. Please add its relative mass with respect to the Sun.")
        relative_masses[planet] = float(input("Relative mass (float number, add the denominator, e.g., 328900.56 for Earth+Moon): \n"))
    
    t_r = float(input("Reference time (JD) = "))
    x = float(input("x = "))
    y = float(input("y = "))
    z = float(input("z = "))
    x_speed = float(input("x_speed = "))
    y_speed = float(input("y_speed = "))
    z_speed = float(input("z_speed = "))

    return x, y, z, x_speed, y_speed, z_speed, planet, t_r

#                   --- Main function (CLI) ---

def main():
    parser = argparse.ArgumentParser(description="Calculate orbital elements of a celestial object.")
    parser.add_argument("--object", type=str, help="Name of the planet (e.g., 'Earth+Moon').")
    parser.add_argument("--t_r", type=float, help="Reference time in Julian days (JD).")
    parser.add_argument("--x", type=float, help="x coordinate in au.")
    parser.add_argument("--y", type=float, help="y coordinate in au.")
    parser.add_argument("--z", type=float, help="z coordinate in au.")
    parser.add_argument("--x_speed", type=float, help="x speed in au/dia.")
    parser.add_argument("--y_speed", type=float, help="y speed in au/dia.")
    parser.add_argument("--z_speed", type=float, help="z speed in au/dia.")
    args = parser.parse_args()

    print("\t########## ORBITAL ELEMENTS CALCULATION ##########\n")

    #   --- Initialize the orbital elements ---
    if args.object is None:
        x, y, z, x_speed, y_speed, z_speed, planet, t_r = manual_input()

    else:
        object = args.object
        t_r = args.t_r
        x = args.x
        y = args.y
        z = args.z
        x_speed = args.x_speed
        y_speed = args.y_speed
        z_speed = args.z_speed

    orbital_objects = OrbitalObjects()

    position, r_magnitude = orbital_objects.position_vector(x, y, z)
    velocity, v_magnitude = orbital_objects.velocity_vector(x_speed, y_speed, z_speed)
    angular_momentum, h_magnitude = orbital_objects.angular_momentum_vector(position, velocity)
    e_vector, e_magnitude = orbital_objects.eccentricity_vector(planet, position, velocity, r_magnitude, angular_momentum)
    a = orbital_objects.semi_major_axis(planet, v_magnitude, r_magnitude)
    dr_dt = orbital_objects.rate_of_change_of_position(r_magnitude, position, velocity) # Auxiliary variable 
    M = orbital_objects.mean_anomaly(a, e_magnitude, r_magnitude, dr_dt)
    i = orbital_objects.inclination(h_magnitude, angular_momentum)
    Omega = orbital_objects.longitude_of_ascending_node(angular_momentum)
    omega = orbital_objects.argument_of_periapsis(e_vector, angular_momentum, h_magnitude)
    t_0 = orbital_objects.proper_time_of_periapsis_passage


if __name__ == "__main__":
    main()