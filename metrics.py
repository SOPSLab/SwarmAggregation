# Project:     SwarmAggregation
# Filename:    metrics.py
# Authors:     Joshua J. Daymude (jdaymude@asu.edu) and Noble C. Harasha
#              (nharasha1202@gmail.com).

"""
metrics: A library of aggregation metrics.
"""

from math import hypot
import numpy as np
from scipy.spatial import ConvexHull, distance_matrix
from welzl import welzl


def sed_circumference(config):
    """
    Takes as input an N x 3 array of robot position and orientation data and
    returns the circumference of the system's smallest enclosing disc.
    """
    _, _, radius = welzl(config[:,:2])
    return 2*np.pi * radius


def hull_perimeter(config):
    """
    Takes as input an N x 3 array of robot position and orientation data and
    returns the perimeter of the system's convex hull.
    """
    hull = ConvexHull(config[:,:2])
    perimeter = 0
    for i in range(len(hull.vertices)):
        v1, v2 = config[hull.vertices[i-1]][:2], config[hull.vertices[i]][:2]
        perimeter += hypot(*(v1 - v2))

    return perimeter


def dispersion(config):
    """
    Takes as input an N x 3 array of robot position and orientation data and
    returns the system's dispersion.
    """
    xs, ys = config[:,0], config[:,1]
    return np.sum(np.sqrt((xs - np.mean(xs))**2 + (ys - np.mean(ys))**2))


def cluster_fraction(config, r, eps=0.05):
    """
    Takes as input an N x 3 array of robot position and orientation data and the
    radius of each robot and returns the fraction of robots in the system's
    largest connected cluster, where 'connected' is within epsilon% touching.
    """
    def DFS(config, i, r, dists, visited, cluster):
        visited[i] = 1
        cluster.append(i)
        for j in range(len(config)):
            if visited[j] == 0 and dists[i][j] <= (2 + eps)*r:
                DFS(config, j, r, dists, visited, cluster)

    dists = distance_matrix(config[:,:2], config[:,:2])
    visited = np.zeros(len(config))
    clusters = []

    for i in range(len(config)):
        if visited[i] == 0:
            cluster = []
            DFS(config, i, r, dists, visited, cluster)
            clusters.append(cluster)

    return max([len(cluster) for cluster in clusters]) / len(config)
