"""
Recomendation: 
    Install VPython for better vizualization 
        pip install vpython
        pip install "setuptools<70"
        
    If you can't install VPython, you can use Matplotlib for 3D plotting as shown in the code.

Usage: 
    python iteration.py --object Mars --time 2451545.0 --steps 1000

Description:
    This code visualizes the orbit of a celestial object (e.g., Mars) based on its orbital elements 
    and the time of observation.
    It calculates the position of the object at different time steps and plots its orbit in a 3D space using Matplotlib
    or VPython.

    # Manual input:
    You can add the parameters manually if you don't want to use command-line arguments.

    # Command-line arguments (CLI):
    You can use the command-line arguments
    --object: Name of the celestial object (default: Mars)
    --time: Time of observation in Julian Date (JD) (default: 2451545.0)
    --steps: Number of steps for plotting the orbit (default: 1000)

"""

from two_body_problem import TwoBodyProblem, orbital_elements, relative_masses, reference_time, K
import numpy as np
import matplotlib.pyplot as plt
import argparse

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
        _, _, _, _, x, y, z = object_solution.general_solution(t)
        x_values.append(x)
        y_values.append(y)
        z_values.append(z)

    try: 
        from vpython import sphere, vector, color, rate, canvas

        # Plot the orbit in 3D using VPython
        scene = canvas(title=f'Orbit of {object_name} from JD {time} to JD {time + T:.2f}', background=color.black)
        sun = sphere(pos=vector(0,0,0), radius=0.1, color=color.yellow, emissive=True)

        object_3d = sphere(radius=0.06, color=color.red, make_trail=True)

        while True:
            for x, y, z in zip(x_values, y_values, z_values):
                object_3d.pos = vector(x, y, z)
                rate(30)  # For smoother animation use 60 fps or higher

    except ImportError:

        # Plot the orbit usign Matplotlib
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot(x_values, y_values, z_values)
        ax.set_xlabel('X (au)')
        ax.set_ylabel('Y (au)')
        ax.set_zlabel('Z (au)')
        ax.set_title(f'Orbit of {object_name} from JD {time} to JD {time + T:.2f}')
        plt.show()

def manual_input():
    print("\t--- MANUAL INPUT ---")
    object_name = input("Enter the name of the celestial object (e.g., Mars): ")
    time = float(input("Enter the time of observation in Julian Date (JD) (e.g., 2451545.0): "))
    steps = int(input("Enter the number of steps for plotting the orbit (e.g., 1000): "))
    
    return object_name, time, steps

def main():

    #   --- Command-line arguments (CLI) OPTIONAL ---
    parser = argparse.ArgumentParser(description='Plot the orbit of a celestial object.')
    parser.add_argument('--object', type=str, help='Name of the celestial object (default: Mars)')
    parser.add_argument('--time', type=float, help='Time of observation in Julian Date (JD) (default: 2451545.0)')
    parser.add_argument('--steps', type=int, help='Number of steps for plotting the orbit (default: 1000)')
    args = parser.parse_args()

    print("\t########## VISUALIZATION OF THE ORBIT ON THE TWO-BODY PROBLEM ########## \n")

    object_name = args.object
    time = args.time
    steps = args.steps

    if not all([object_name, time, steps]):
        object_name, time, steps = manual_input()

    plot_orbit(object_name, time, steps)

if __name__ == "__main__":
    main()