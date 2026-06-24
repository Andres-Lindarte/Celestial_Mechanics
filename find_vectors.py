"""
Usage:
    python find_vectors.py --planet "Earth+Moon" --x -0.1771354615 --y 0.9672416229 --z -0.0000039000 --x_speed -0.0172020989 --y_speed -0.0034906585 --z_speed 0.0000000000

Description:
    This script calculates important vectors and quantities in celestial mechanics, 
    given the position and velocity vectors of the planet.
    The vectors and quantities calculated include:
    - Position vector
    - Velocity vector
    - Angular momentum vector
    - Cenital angle
    - Rate of change of position
    - Rate of change of true anomaly
    - Eccentricity vector
    - Total energy per unit mass

    # Manual input:
    You can give the position and velocity coordinates manually by running the script without arguments.

    # Command-line arguments:
    You can use the command-line arguments
    --planet: Name of the planet (e.g., 'Earth+Moon').
    --x: x coordinate in ua.
    --y: y coordinate in ua.
    --z: z coordinate in ua.
    --x_speed: x speed in ua/dia.
    --y_speed: y speed in ua/dia.
    --z_speed: z speed in ua/dia.

"""


from useful import UsefulFunctions, K, RELATIVE_MASSES 
import numpy as np
import argparse

#                    --- Useful functions ---
    
useful = UsefulFunctions()

#                    --- Finding the cuantities ---

class FindVectors:

    """
    This class contains methods to calculate important vectors and quantities in celestial mechanics,
    such as:
    - Position vector
    - Velocity vector
    - Angular momentum vector
    - Cenital angle
    - Rate of change of position
    - Rate of change of true anomaly
    - Eccentricity vector
    - Total energy per unit mass
    
    """

    @staticmethod
    def position_vector (x,y,z):
        # Calculate the magnitude of the position vector
        position = np.array([x, y, z])
        r_magnitude = np.linalg.norm(position)
        print(f"Position: (->r)=[{x}, {y}, {z}]     (r)= {r_magnitude:.8f} ua \n")
        
        return position, r_magnitude

    @staticmethod
    def velocity_vector (x_speed, y_speed, z_speed):
        # Calculate the magnitude of the velocity vector
        velocity = np.array([x_speed, y_speed, z_speed])
        v_magnitude = np.linalg.norm(velocity)
        print(f"Velocity: (->v)=[{x_speed}, {y_speed}, {z_speed}]     (v)= {v_magnitude:.8f} ua/dia \n")
        
        return velocity, v_magnitude

    @staticmethod
    def angular_momentum_vector (position, velocity):
        # Calculate the cross product of position and velocity vectors (the angular momentum vector per unit mass h)
        angular_momentum = np.cross(position, velocity)
        h_magnitude = np.linalg.norm(angular_momentum)
        print(f"Angular Momentum per unit mass: (->h)=[{angular_momentum[0]:.8f}, {angular_momentum[1]:.8f}, {angular_momentum[2]:.8f}]     \n (h)= {h_magnitude:.8f} ua^2/dia \n")

        return angular_momentum, h_magnitude

    @staticmethod
    def cenital_angle (r_magnitude, v_magnitude, position, velocity):
        # Calculate the cenital angle (cen) using the formula: cen = arccos((position . velocity )/(rv))
        cen = np.arccos((np.dot(position, velocity)) / (r_magnitude * v_magnitude))
        sign, d, m, s = useful.radians_to_dms(cen)
        print(f"Cenital Angle: (cen)= {np.degrees(cen):.8f} = {sign}{d}° {m}' {s:.2f}\" \n")
        return cen

    @staticmethod
    def rate_of_change_of_position (r_magnitude, position, velocity):
        # Calculate the rate of change of position (dr/dt) using the formula: dr/dt = (position . velocity) / r
        dr_dt = np.dot(position, velocity) / r_magnitude
        print(f"Rate of Change of Position: (dr/dt)= {dr_dt:.8f} ua/dia \n")
        return dr_dt

    @staticmethod
    def rate_of_change_of_true_anomaly (r_magnitude, h_magnitude,):
        # Calculate the rate of change of true anomaly (dtheta/dt) using the formula: dtheta/dt = (h) / (r^2)
        dtheta_dt = h_magnitude / (r_magnitude**2)
        print(f"Rate of Change of True Anomaly: (dtheta/dt)= {dtheta_dt:.8f} rad/dia \n")
        return dtheta_dt

    @staticmethod
    def eccentricity_vector (planet, position, velocity, r_magnitude, angular_momentum):
        # Calculate the eccentricity vector (e) using the formula: e = -(h x v) / mu - (->r) / r
        mu = K**2 * (1 + RELATIVE_MASSES[planet])  # Standard gravitational parameter for the planet
        e_vector = -(np.cross(angular_momentum, velocity) / mu) - (position / r_magnitude)
        e_magnitude = np.linalg.norm(e_vector)
        print(f"Eccentricity Vector: (->e)=[{e_vector[0]:.8f}, {e_vector[1]:.8f}, {e_vector[2]:.8f}]     \n (e)= {e_magnitude:.8f} \n")
        return e_vector, e_magnitude

    @staticmethod
    def total_energy_per_unit_mass (planet, v_magnitude, r_magnitude):
        # Calculate the total energy per unit mass (E) using the formula: E_mass = (v^2)/2 - mu/r
        mu = K**2 * (1 + RELATIVE_MASSES[planet])  # Standard gravitational parameter for the planet
        E_mass = (v_magnitude**2) / 2 - mu / r_magnitude
        print(f"Total Energy per Unit Mass: (E)= {E_mass:.8f} ua^2/dia^2 \n")
        return E_mass

#                    --- Useful functions (Manual input) ---

def manual_input():

    print("\t--- MANUAL INPUT ---")
    print("Choose the planet for the calculations:")
    planet = input("Planet: ")
    if planet not in RELATIVE_MASSES:
        print("Invalid planet name. Please choose from the following: Mercury, Venus, Earth+Moon, Mars, Jupiter, Saturn, Uranus, Neptune.")
        exit()

    print("Enter the position vector (x, y, z) in astronomical units (ua):")
    x = float(input("x = "))
    y = float(input("y = "))
    z = float(input("z = "))

    print("\nEnter the velocity vector (x, y, z) in astronomical units (ua/dia):")
    x_speed = float(input("x_speed = "))
    y_speed = float(input("y_speed = "))
    z_speed = float(input("z_speed = "))

    return planet, x, y, z, x_speed, y_speed, z_speed

#                   --- Main function ---

def main():

    parser = argparse.ArgumentParser(description="Calculate important vectors and quantities in celestial mechanics.")
    parser.add_argument('--planet', type=str, help='Name of the planet (e.g., Mercury, Venus, Earth+Moon, Mars, Jupiter, Saturn, Uranus, Neptune)')
    parser.add_argument('--x', type=float, help='x coordinate in ua')
    parser.add_argument('--y', type=float, help='y coordinate in ua')
    parser.add_argument('--z', type=float, help='z coordinate in ua')
    parser.add_argument('--x_speed', type=float, help='x speed in ua/dia')
    parser.add_argument('--y_speed', type=float, help='y speed in ua/dia')
    parser.add_argument('--z_speed', type=float, help='z speed in ua/dia')
    args = parser.parse_args()

    if args.planet not in RELATIVE_MASSES:
        print("Invalid planet name. Please choose from the following: Mercury, Venus, Earth+Moon, Mars, Jupiter, Saturn, Uranus, Neptune.")
        exit()

    print("\t########## IMPORTANT VECTORS AND QUANTITIES ##########\n")
    
    #   --- Initialize the requirements ---
    if args.planet is None:
        planet, x, y, z, x_speed, y_speed, z_speed = manual_input()
    else:
        planet = args.planet
        x = args.x
        y = args.y
        z = args.z
        x_speed = args.x_speed
        y_speed = args.y_speed
        z_speed = args.z_speed

    #   --- Calculate the vectors and quantities ---
    print("\n\t--- RESULTS ---\n")
    
    find_vectors = FindVectors()

    position, r_magnitude = find_vectors.position_vector(x, y, z)
    velocity, v_magnitude = find_vectors.velocity_vector(x_speed, y_speed, z_speed)
    angular_momentum, h_magnitude = find_vectors.angular_momentum_vector(position, velocity)
    theta = find_vectors.cenital_angle(r_magnitude, v_magnitude, position, velocity)
    dr_dt = find_vectors.rate_of_change_of_position(r_magnitude, position, velocity)
    dtheta_dt = find_vectors.rate_of_change_of_true_anomaly(r_magnitude, h_magnitude)   
    e_vector, e_magnitude = find_vectors.eccentricity_vector(planet, position, velocity, r_magnitude, angular_momentum)
    E_mass = find_vectors.total_energy_per_unit_mass(planet, v_magnitude, r_magnitude)

if __name__ == "__main__":
    main()