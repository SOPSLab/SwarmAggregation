# Project:     SwarmAggregation
# Filename:    aggregation.py
# Authors:     Joshua J. Daymude (jdaymude@asu.edu) and Noble C. Harasha
#              (nharasha1202@gmail.com).

"""
aggregation: A physical simulation of e-puck robots performing aggregation
             according to the Gauci et al. (2014) algorithm.
"""

from metrics import dispersion
import numpy as np
from tqdm import tqdm


def init_rand(N, r, rng):
    """
    Initialize N robots with radii r uniformly at random in a box centered at
    (0,0) such that the robots occupy roughly 2% of the space.
    """
    config = np.zeros((N, 3))
    arena_side_len = np.sqrt(N * np.pi * r**2 / 0.02)

    for i in range(N):
        while True:
            # Propose a random position for the i-th robot and accept it if it
            # does not overlap with any existing ones.
            config[i][0] = (rng.random() - 0.5) * arena_side_len  # X
            config[i][1] = (rng.random() - 0.5) * arena_side_len  # Y
            config[i][2] = rng.random() * 2*np.pi                 # Theta
            if np.all(np.linalg.norm(config[:i,:2]-config[i][:2], axis=1) > 2*r):
                break

    return config


def init_symm(N, R):
    """
    Initialize N robots with rotational radii R in a symmetric cycle.
    """
    config = np.zeros((N, 3))
    cycle_radius = (N + 1) * R
    cycle_step = 2*np.pi / N

    for i in range(N):
        config[i][0] = cycle_radius * np.cos(cycle_step * i)  # X
        config[i][1] = cycle_radius * np.sin(cycle_step * i)  # Y
        config[i][2] = cycle_step * i                         # Theta

    return config


def ideal(N, r):
    """
    Takes as input a number of robots N with radii r and returns the hexagonal
    packing configuration.
    """
    # First compute the robot positions in a hexagon in triangular coordinates.
    points_tri, layer = [[0,0]], 1
    while len(points_tri) < N:
        # Start a new layer at (layer, 0).
        x, y = layer, 0
        points_tri.append([x,y])
        if len(points_tri) == N:
            break

        # Extend along the six segments of the layer.
        while y > -layer and len(points_tri) < N:  # Down-left.
            y -= 1
            points_tri.append([x,y])
        while x > 0 and len(points_tri) < N:       # Left.
            x -= 1
            points_tri.append([x,y])
        while y < 0 and len(points_tri) < N:       # Up-left.
            x -= 1
            y += 1
            points_tri.append([x,y])
        while y < layer and len(points_tri) < N:   # Up-right.
            y += 1
            points_tri.append([x,y])
        while x < 0 and len(points_tri) < N:       # Right.
            x += 1
            points_tri.append([x,y])
        while y > 1 and len(points_tri) < N:       # Down-right.
            x += 1
            y -= 1
            points_tri.append([x,y])
        layer += 1

    # Convert triangular coordinates into Cartesian ones.
    config = np.zeros((N, 3))
    for i in range(N):
        x_tri, y_tri = points_tri[i]
        config[i][:2] = [r * (x_tri + (y_tri / 2)), r * (np.sqrt(3)/2) * y_tri]

    return config


def spring_constant(R, r, m, w0, step):
    """
    Calculate the spring constant for collisions such that the resulting spring
    force balances the maximum possible overlap.

    Inputs:
    - R (float): distance from a robot's center to its center of rotation (m)
    - r (float): radius of a robot (m)
    - m (float): mass of a robot (kg)
    - w0 (float): rot. speed of a robot about its center of rotation (rad/s)
    - step (float): wall-clock duration of a time step (s)

    Returns: K (double) the spring constant
    """
    min_dist = np.sqrt(8 * R**2 * (1 - np.cos(w0 * step)) + 4 * r**2 \
                       + 8 * R * r * np.sin(w0 * step))
    return (min_dist * m) / ((2 * r - min_dist) * step**2)


def sense(config, i, r, sensor):
    """
    Use the line/cone-of-sight sensor to look for other robots.

    Inputs:
    - config: N x 3 array of robot position and orientation data
    - i (int): index of the robot doing the sensing
    - r (float): radius of a robot (m)
    - sensor (float): size of the line/cone-of-sight sensor (rad)

    Returns: True if another robot is in the sensor region; False otherwise.
    """
    x_i, y_i, theta_i = config[i]
    cw, ccw = theta_i - sensor / 2, theta_i + sensor / 2

    for j in np.arange(len(config)):
        x_j, y_j, _ = config[j]

        # Remove the half-plane that is "more clockwise" than the "more
        # clockwise" line that is parallel to the clockwise boundary and tangent
        # to robot i.
        if cw < np.pi/2 or cw > 3*np.pi/2:
            if y_j < np.tan(cw) * (x_j - x_i) + y_i - r / np.cos(cw):
                continue
        elif cw > np.pi/2 and cw < 3*np.pi/2:
            if y_j > np.tan(cw) * (x_j - x_i) + y_i - r / np.cos(cw):
                continue
        elif np.isclose(cw, np.pi/2):
            if x_j > x_i + r:
                continue
        else:  # np.isclose(cw, 3*np.pi/2)
            if x_j < x_i - r:
                continue

        # Remove the half-plane that is "more counter-clockwise" than the "more
        # counter-clockwise" line that is parallel to the clockwise boundary and
        # tangent to robot i.
        if ccw < np.pi/2 or ccw > 3*np.pi/2:
            if y_j > np.tan(ccw) * (x_j - x_i) + y_i + r / np.cos(ccw):
                continue
        elif ccw > np.pi/2 and ccw < 3*np.pi/2:
            if y_j < np.tan(ccw) * (x_j - x_i) + y_i + r / np.cos(ccw):
                continue
        elif np.isclose(ccw, np.pi/2):
            if x_j < x_i - r:
                continue
        else:  # np.isclose(ccw, 3*np.pi/2)
            if x_j > x_i + r:
                continue

        # Remove the half-plane "behind" robot i (perpendicular to the sensor).
        # BUG: This both excludes things that are in fact seen and includes
        # things that are in fact not seen.
        if theta_i > 0 and theta_i < np.pi:
            if y_j <= np.tan(theta_i + np.pi/2) * (x_j - x_i) + y_i:
                continue
        elif theta_i > np.pi and theta_i < 2*np.pi:
            if y_j >= np.tan(theta_i + np.pi/2) * (x_j - x_i) + y_i:
                continue
        elif np.isclose(theta_i, 0) or np.isclose(theta_i, 2*np.pi):
            if x_j <= x_i:
                continue
        else:  # np.isclose(theta_i, np.pi)
            if x_j >= x_i:
                continue

        # If none of the above occur, then robot j is seen by robot i.
        return True

    return False


def update(config, R, r, m, w0, w1, K, sensor, noise, step, rng):
    """
    Compute robot positions and orientations after a single time step.

    Inputs:
    - config: N x 3 array of robot position and orientation data
    - R (float): distance from a robot's center to its center of rotation (m)
    - r (float): radius of a robot (m)
    - m (float): mass of a robot (kg)
    - w0 (float): rot. speed of a robot about its center of rotation (rad/s)
    - w1 (float): rot. speed of a robot in place (rad/s)
    - K (float): spring constant for collisions (N/m)
    - sensor (float): size of the line/cone-of-sight sensor (rad)
    - noise ((str, float)): ('err', p) for error probability with probability p
                            ('mot', f) for motion noise with force f (N)
    - step (float): wall-clock duration of a time step (s)
    - rng: random number generator

    Returns: updated N x 3 array of robot position and orientation data
    """
    next = np.copy(config)

    # Compute collision forces.
    forces = np.zeros((len(config), 2))
    for i in range(len(config)):
        for j in range(i+1, len(config)):
            # If robots i and j overlap, then apply a spring force to both.
            delta = config[i][:2] - config[j][:2]
            dist = np.linalg.norm(delta)
            if dist <= 2*r:
                forces[i] += (delta / dist) * K * (2*r - dist)
                forces[j] -= (delta / dist) * K * (2*r - dist)

    # Compute motion noise forces.
    if noise[0] == 'mot':
        for i in range(len(config)):
            force = rng.random() * noise[1]
            theta = rng.random() * 2*np.pi
            forces[i] += force * np.array([np.cos(theta), np.sin(theta)])

    # Compute translation and rotation from algorithm drive.
    for i in range(len(config)):
        # Use line/cone-of-sight and invert response if using error probability.
        see_other = sense(config, i, r, sensor)
        if noise[0] == 'err' and rng.random() < noise[1]:
            see_other = not see_other

        # Compute movement updates based on the sensing response.
        if see_other:  # Rotate in place.
            next[i][2] = (config[i][2] + (w1 * step) + (2*np.pi)) % (2*np.pi)
        else:  # Rotate clockwise around center of rotation.
            # Compute location of the center of rotation, c.
            c_angle = (config[i][2] + np.pi/2)
            c = config[i][:2] + R * np.array([np.cos(c_angle), np.sin(c_angle)])

            # Compute rotation around c.
            c_angle = (c_angle + np.pi + (w0 * step))
            next[i][:2] = c + R * np.array([np.cos(c_angle), np.sin(c_angle)])
            next[i][2] = (config[i][2] + (w0 * step) + (2*np.pi)) % (2*np.pi)

    # Integrate forces and apply to the updated configuration.
    next[:,:2] += forces * step**2 / m

    return next


def aggregation(N=50, R=0.1445, r=0.037, m=0.152, w0=-0.75, w1=-5.02, sensor=0,\
                noise=('err', 0), time=60, step=0.005, stop=None, init='rand', \
                seed=None, silent=False):
    """
    Execute an aggregation experiment.

    Throughout, configuration data tracks:
    - X (float): the x-position of each robot (m)
    - Y (float): the y-position of each robot (m)
    - Theta (float): the orientation of each robot's sight sensor (rad)

    Inputs:
    - N (int): number of robots
    - R (float): distance from a robot's center to its center of rotation (m)
    - r (float): radius of a robot (m)
    - m (float): mass of a robot (kg)
    - w0 (float): rot. speed of a robot about its center of rotation (rad/s)
    - w1 (float): rot. speed of a robot in place (rad/s)
    - sensor (float): size of the line/cone-of-sight sensor (rad)
    - noise ((str, float)): ('err', p) for error probability with probability p
                            ('mot', f) for motion noise with force f (N)
    - time (float): wall-clock duration of experiment (s)
    - step (float): wall-clock duration of a time step (s)
    - stop (float): if not None, experiment will stop before 'time' if the
                    system's dispersion is within stop% of the ideal value
    - init (str): 'rand' for random, 'symm' for symmetric
    - seed (int): random seed
    - silent (bool): False if progress should be shown on command line

    Returns (history, i_f):
    - history: #steps x N x 3 array of robot X/Y/Theta over time
    - i_f (int): index of history containing final configuration
    """
    # Initialize the random number generation.
    rng = np.random.default_rng(seed)

    # Initialize the data array. Each robot has (X, Y, Theta) at each step.
    steps = np.arange(0, time + step, step)
    history = np.zeros((len(steps), N, 3))

    # Initialize the robots according to the specified method.
    if init == 'rand':
        history[0] = init_rand(N, r, rng)
    elif init == 'symm':
        history[0] = init_symm(N, r)
    else:
        assert False, 'ERROR: Unrecogized initialization method: ' + init

    # For efficiency's sake, pre-calculate the ideal dispersion for N robots and
    # the spring constant K for collisions.
    disp_ideal = dispersion(ideal(N, r))
    K = spring_constant(R, r, m, w0, step)

    # Simulate the experiment duration by time step.
    for i, t in enumerate(tqdm(steps[1:], disable=silent)):
        # Compute updates to robot positions and orientations.
        history[i+1] = update(history[i], R, r, m, w0, w1, K, sensor, noise, \
                              step, rng)

        # If specified, check the stopping condition once per second.
        if stop != None and float.is_integer(t) and \
           dispersion(history[i+1]) <= (1 + stop) * disp_ideal:
            return (history, i+1)

    # Return history of configurations.
    return (history, -1)
