from two_body_problem import TwoBodyProblem, orbital_elements, relative_masses, reference_time, K
import numpy as np
import matplotlib.pyplot as plt

def plot_orbit(object_name, time, steps):

    """
    This function plots the orbit of a celestial object based on its orbital elements and the time of observation.
    It uses the TwoBodyProblem class to calculate the position of the object at the given time and then plots
    its orbit in a 3D space.
    
    Parameters:
    - object_name: Name of the celestial object (e.g., "Mars")
    - time: Time of observation in Julian Date (JD)
    - orbital_elements: A dictionary containing the orbital elements of the object, including:
        - M_r: Mean anomaly at reference time
        - a: Semi-major axis
        - e: Eccentricity
        - i: Inclination
        - Omega: Longitude of the ascending node
        - omega: Argument of periapsis
    
    Returns:
    - A 3D plot of the object's orbit.
    
    """
    x_values , y_values, z_values = [], [], []

    # Calculates the period of the object using the formula: T = 2 * pi * sqrt(a^3 / (K^2 * (1 + relative_mass)))
    T = 2 * np.pi * np.sqrt(orbital_elements[object_name]["a"]**3 / (K**2 * (1 + relative_masses[object_name])))

    object_solution = TwoBodyProblem(
        object_name=object_name, 
        orbital_elements=orbital_elements[object_name], 
        relative_mass=relative_masses.get(object_name, 0), 
        reference_time=reference_time,
        verbose=False)
    
    for t in np.linspace(time, time + T, steps):
        # Calculate the position of the object at time t
        M, E, r, theta, x, y, z = object_solution.general_solution(t)
        x_values.append(x)
        y_values.append(y)
        z_values.append(z)

    # Plot the orbit in 3D
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(x_values, y_values, z_values)
    ax.set_xlabel('X (au)')
    ax.set_ylabel('Y (au)')
    ax.set_zlabel('Z (au)')
    ax.set_title(f'Orbit of {object_name} from JD {time} to JD {time + T:.2f}')
    plt.show()

def main():
    object_name = "Mars"
    time = 2451545.0  # Example Julian Date (JD)
    steps = 1000  # Number of steps for plotting the orbit

    plot_orbit(object_name, time, steps)

if __name__ == "__main__":
    main()