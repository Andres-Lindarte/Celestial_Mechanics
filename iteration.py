from two_body_problem import TwoBodyProblem
import numpy as np


def plot_orbit(object_name, time, orbital_elements):

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
    pass

