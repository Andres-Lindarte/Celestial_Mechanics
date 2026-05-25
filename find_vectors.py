import numpy as np

#                    --- Useful funtions ---
def radians_to_dms(rad):
    deg = np.degrees(rad)
    sign_string = "-" if deg < 0 else "+"
    deg = abs(deg)
    d = int(deg)
    m = int((deg - d) * 60)
    s = (deg - d - m/60) * 3600
    return sign_string, d, m, s
    
def position_vector (x,y,z):
    # Calculate the magnitude of the position vector
    position = np.array([x, y, z])
    r_magnitude = np.linalg.norm(position)
    print(f"Position: (->r)=[{x}, {y}, {z}]     (r)= {r_magnitude:.8f} ua \n")
    
    return position, r_magnitude

def velocity_vector (x_speed, y_speed, z_speed):
    # Calculate the magnitude of the velocity vector
    velocity = np.array([x_speed, y_speed, z_speed])
    v_magnitude = np.linalg.norm(velocity)
    print(f"Velocity: (->v)=[{x_speed}, {y_speed}, {z_speed}]     (v)= {v_magnitude:.8f} ua/dia \n")
    
    return velocity, v_magnitude

def angular_momentum_vector (position, velocity):
    # Calculate the cross product of position and velocity vectors (the angular momentum vector per unit mass h)
    angular_momentum = np.cross(position, velocity)
    h_magnitude = np.linalg.norm(angular_momentum)
    print(f"Angular Momentum per unit mass: (->h)=[{angular_momentum[0]:.8f}, {angular_momentum[1]:.8f}, {angular_momentum[2]:.8f}]     \n (h)= {h_magnitude:.8f} ua^2/dia \n")

    return angular_momentum, h_magnitude

def cenital_angle (r_magnitude, v_magnitude, position, velocity):
    # Calculate the cenital angle (cen) using the formula: cen = arccos((position . velocity )/(rv))
    cen = np.arccos((np.dot(position, velocity)) / (r_magnitude * v_magnitude))
    sign, d, m, s = radians_to_dms(cen)
    print(f"Cenital Angle: (cen)= {np.degrees(cen):.8f} = {sign}{d}° {m}' {s:.2f}\" \n")
    return cen

def rate_of_change_of_position (r_magnitude, position, velocity):
    # Calculate the rate of change of position (dr/dt) using the formula: dr/dt = (position . velocity) / r
    dr_dt = np.dot(position, velocity) / r_magnitude
    print(f"Rate of Change of Position: (dr/dt)= {dr_dt:.8f} ua/dia \n")
    return dr_dt

def rate_of_change_of_true_anomaly (r_magnitude, h_magnitude,):
    # Calculate the rate of change of true anomaly (dtheta/dt) using the formula: dtheta/dt = (h) / (r^2)
    dtheta_dt = h_magnitude / (r_magnitude**2)
    print(f"Rate of Change of True Anomaly: (dtheta/dt)= {dtheta_dt:.8f} rad/dia \n")
    return dtheta_dt

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

k = 0.01720209908  # Gaussian gravitational constant in au^(3/2)/dia

def eccentricity_vector (planet, position, velocity, r_magnitude, angular_momentum):
    # Calculate the eccentricity vector (e) using the formula: e = -(h x v) / mu - (->r) / r
    mu = k**2 * (1 + relative_masses[planet])  # Standard gravitational parameter for the planet
    e_vector = -(np.cross(angular_momentum, velocity) / mu) - (position / r_magnitude)
    e_magnitude = np.linalg.norm(e_vector)
    print(f"Eccentricity Vector: (->e)=[{e_vector[0]:.8f}, {e_vector[1]:.8f}, {e_vector[2]:.8f}]     \n (e)= {e_magnitude:.8f} \n")
    return e_vector, e_magnitude

def total_energy_per_unit_mass (planet, v_magnitude, r_magnitude):
    # Calculate the total energy per unit mass (E) using the formula: E_mass = (v^2)/2 - mu/r
    mu = k**2 * (1 + relative_masses[planet])  # Standard gravitational parameter for the planet
    E_mass = (v_magnitude**2) / 2 - mu / r_magnitude
    print(f"Total Energy per Unit Mass: (E)= {E_mass:.8f} ua^2/dia^2 \n")
    return E_mass

#                   --- Main function ---
print("\t--- IMPORTANT VECTORS ---\n")
print("Choose the planet for the calculations:")
planet = input("Planet: ")
if planet not in relative_masses:
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

print("\n\t--- RESULTS ---\n")
position, r_magnitude = position_vector(x, y, z)
velocity, v_magnitude = velocity_vector(x_speed, y_speed, z_speed)
angular_momentum, h_magnitude = angular_momentum_vector(position, velocity)
theta = cenital_angle(r_magnitude, v_magnitude, position, velocity)
dr_dt = rate_of_change_of_position(r_magnitude, position, velocity)
dtheta_dt = rate_of_change_of_true_anomaly(r_magnitude, h_magnitude)   
e_vector, e_magnitude = eccentricity_vector(planet, position, velocity, r_magnitude, angular_momentum)
E_mass = total_energy_per_unit_mass(planet, v_magnitude, r_magnitude)