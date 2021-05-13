# Project:     SwarmAggregation
# Filename:    aggregation.py
# Authors:     Joshua J. Daymude (jdaymude@asu.edu) and Noble C. Harasha
#              (nharasha1202@gmail.com).

"""
metrics: A library of aggregation metrics.
"""

import numpy as np


def sed_circumference(config):
    """
    Takes as input an N x 3 array of robot position and orientation data and
    returns the circumference of the system's smallest enclosing disc.
    """
    pass  # TODO.


def convex_hull(config):
    """
    Takes as input an N x 3 array of robot position and orientation data and
    returns the perimeter of the system's convex hull.
    """
    pass  # TODO.


def dispersion(config):
    """
    Takes as input an N x 3 array of robot position and orientation data and
    returns the system's dispersion.
    """
    xs, ys = config[:,0], config[:,1]
    return np.sum(np.sqrt((xs - np.mean(xs))**2 + (ys - np.mean(ys))**2))    


def cluster_fraction(config):
    """
    Takes as input an N x 3 array of robot position and orientation data and
    returns the system's cluster fraction.
    """
    pass  # TODO.
