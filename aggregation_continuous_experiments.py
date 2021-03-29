# Project:     Swarm Aggregation
# Filename:    aggregation_continuous_experiments.py
# Authors:     Noble C. Harasha (nharasha1202@gmail.com)

"""
experiments: A framework for defining and running experiments using the swarm
             aggregation simulations.
"""

import argparse
from itertools import product
import pickle
from  aggregation_continuous import *
import datetime
import math
import matplotlib.animation
import matplotlib.cm as cm
from matplotlib.collections import LineCollection
from matplotlib.collections import PatchCollection
from matplotlib.collections import PolyCollection
from matplotlib.colors import LogNorm
from matplotlib.colors import to_hex
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
import random
import statistics as stats       # for statistics.mean(), .pstdev()
















class Circle:
    def __init__(self, C, R):
        self.C = C
        self.R = R

def isInsideCircle(circle, point):
    return distArrs(circle.C, point) <= circle.R

def circleFromThree(a, b, c):
    mpAB = [ (a[0] + b[0]) / 2.0, (a[1] + b[1]) / 2.0 ]
    mpBC = [ (b[0] + c[0]) / 2.0, (b[1] + c[1]) / 2.0 ]
    if (a[0] == b[0] and b[1] == c[1]):
        centerX = mpBC[0]
        centerY = mpAB[1]
    elif (b[0] == c[0] and a[1] == b[1]):
        centerX = mpAB[0]
        centerY = mpBC[1]
    elif (a[0] == b[0]):
        centerY = mpAB[1]
        slopeBC = (c[1] - b[1]) / (c[0] - b[0])
        L2slope = (-1) * (1/slopeBC)
        centerX = ( (centerY - mpBC[1]) / L2slope ) + mpBC[0]
    elif (b[0] == c[0]):
        centerY = mpBC[1]
        slopeAB = (b[1] - a[1]) / (b[0] - a[0])
        L1slope = (-1) * (1/slopeAB)
        centerX = ( (centerY - mpAB[1]) / L1slope) + mpAB[0]
    elif (a[1] == b[1]):
        centerX = mpAB[0]
        slopeBC = (c[1] - b[1]) / (c[0] - b[0])
        L2slope = (-1) * (1/slopeBC)
        centerY = ( L2slope * (centerX - mpBC[0]) ) + mpBC[1]
    elif (b[1] == c[1]):
        centerX = mpBC[0]
        slopeAB = (b[1] - a[1]) / (b[0] - a[0])
        L1slope = (-1) * (1/slopeAB)
        centerY = ( L1slope * (centerX - mpAB[0]) ) + mpAB[1]
    else:
        slopeAB = (b[1] - a[1]) / (b[0] - a[0])
        slopeBC = (c[1] - b[1]) / (c[0] - b[0])
        L1slope = (-1) * (1/slopeAB)
        L2slope = (-1) * (1/slopeBC)
        centerX = ( mpBC[1] - mpAB[1] + (L1slope * mpAB[0]) - (L2slope * mpBC[0]) ) / (L1slope - L2slope)
        centerY = ( L1slope * (centerX - mpAB[0]) ) + mpAB[1]
    center = [centerX, centerY]
    return Circle(center, distArrs(center, b))

def circleFromTwo(a, b):
    center = [ (a[0] + b[0]) / 2.0, (a[1] + b[1]) / 2.0 ]
    return Circle(center, distArrs(a, b) / 2.0)

def isValidCircle(circle, points):
    n = len(points)
    for i in range(n):
        if (isInsideCircle(circle, points[i]) == False):
            return False
            break
    return True

def minCircleTrivial(points):
    l = len(points)
    if (l == 0):
        return Circle([0, 0], 0)
    elif (l == 1):
        return Circle(points[0], 0)
    elif (l == 2):
        return circleFromTwo(points[0], points[1])
    for i in range(3):
        for j in range(i+1, 3):
            c = circleFromTwo(points[i], points[j])
            if (isValidCircle(c, points) == True):
                return c
    return circleFromThree(points[0], points[1], points[2])

def swapListElemPositions(list, pos1, pos2):
    list[pos1], list[pos2] = list[pos2], list[pos1]
    return list

def welzlHelper(points, boundaryPoints, numLeftToProcess):
    if (not boundaryPoints == False):
        if (numLeftToProcess == 0 or len(boundaryPoints) == 3):
            return minCircleTrivial(boundaryPoints)
    idx = random.randint(0, numLeftToProcess-1)
    p = points[idx]
    points = swapListElemPositions(points, idx, numLeftToProcess - 1)
    d = welzlHelper(points, boundaryPoints, numLeftToProcess - 1)
    if (isInsideCircle(d, p) == True):
        return d
    boundaryPoints.append(p)
    return welzlHelper(points, boundaryPoints, numLeftToProcess - 1)

def welzl(points):
    P_copy = points
    random.shuffle(P_copy)
    return welzlHelper(P_copy, list(), len(P_copy))






def DFS(robot, system, clusterArr):
    robot[2] = 1
    clusterArr.append(robot)

    for s in range(len(system)):
        if (distArrs(robot, system[s]) <= 2*3.7 + 0.5 and system[s][2] == 0):
            DFS(system[s], system, clusterArr)










def distArrs(a, b):
    return math.sqrt((b[0] - a[0])**2 + (b[1]-a[1])**2)








class Experiment(object):
    """
    A flexible, unifying framework for experiments.
    """

    def __init__(self, _id, _params={}, _init='random', _noise='errorprob', \
                 _stopping=True, _savehistory=True, _iters=1, _seed=None):
        """
        Takes as input:
        - _id (str): identifier for the experiment, e.g., '1' or 'baseline'.
        - _params (dict): the full parameter set for the aggregation runs.
            {
              'N' (system size) : [int: >= 1],
              'S' (# seconds) : [int: >= 1], *only taken into account if _stopping == 'True'
              'T' (size of cone-of-sight sensor): [float: >= 0 and <= pi],
              'noise_p' (amount of noise) : [float: errorprob: >= 0 and <= 1; motion: >= 0 and <= max],
              'time_step' (# ms for each step; time interval) : [float: > 0] *must sum with itself some integer number of times to 1
            }
        - _init (str): the initialization method for the agents' positions; 'random', 'symmetric', or 'custom'.
        - _noise (str): the noise form; 'errorprob' or 'motion'.
        - _stopping (bool): whether to use a stopping condition or simply run through a specified number of steps; 'True' uses a stopping condition.
        - _savehistory (bool): determines if the position and orientation data of all runs over all steps will be saved; 'False' results in # steps of each run being only saved data
        - _iters (int): the number of iterated runs for each parameter setting.
        - _seed (int): the seed for random number generation.
        """
        # Unpack singular parameters.
        self.id, self.init, self.noise, self.stopping, self.savehistory, self.iters, self.seed = _id, _init, _noise, _stopping, _savehistory, _iters, _seed

        # Unpack aggregation parameters.
        defaults = {'N' : [50], 'S' : [60], 'T' : [0], 'noise_p' : [.05], 'time_step' : [.5]}
        params = [_params[p] if p in _params else defaults[p] for p in defaults]
        self.params = list(product(*params))

        # Set up data and results filenames.
        self.fname = 'exp{}_{}'.format(self.id, self.seed)

        time_of_execution = str(datetime.datetime.now().time())
        punc = '''!()-[]{};:'"\, <>./?@#$%^&*_~'''
        for ele in time_of_execution:
            if ele in punc:
                time_of_execution = time_of_execution.replace(ele, "")

        self.fname = self.fname + '_' + time_of_execution

        # Instantiate a list to hold runs data. This data will have shape
        # R x I x S x N x 3. R is the number of runs (i.e., unique parameter
        # combinations); I is the number of iterations per run; S and N are as
        # they are in the aggregation framework; 3 represents the 3 quantitative
        # attributes/values that are reported for each robot (x position, y
        # position, theta orientation).
        self.runs_data = [[] for p in self.params]



    def run(self):
        tqdm.write('Running Experiment ' + self.id + '...')

        # Set up random seeds for iterated runs.
        rng = np.random.default_rng(self.seed)
        run_seeds = rng.integers(0, 2**32, size=self.iters)

        # For each parameter combination, do iterated runs of aggregation.
        for i, param in enumerate(tqdm(self.params, desc='Simulating runs')):
            # run_seeds = rng.integers(0, 2**32, size=self.iters)
            N, S, T, noise_p, time_step = param
            for seed in tqdm(run_seeds, desc='Iterating run'):
                run_data = aggregation(_N=N, _T=T, _init=self.init, \
                                       _noise=(self.noise, noise_p), \
                                       _time_step=time_step, \
                                       _stopping=(self.stopping, int(S*(1000/time_step))), \
                                       _savehistory=self.savehistory, \
                                       _seed=seed)
                self.runs_data[i].append(run_data)



    def save(self):
        """
        Saves this experiment, including all parameters and run data, to a file
        named according to the experiment's ID and seed.
        """
        tqdm.write('Saving Experiment ' + self.id + '...')
        with open('data/' + self.fname + '.pkl', 'wb') as f:
            pickle.dump(self, f)






    def dispersion(self, _run_configs):
        """
        Takes as input an S x N x 3 run of position/orientation data and returns
        a 1 x S array of the value of the 2nd moment dispersion metric at each
        step.
        """
        S = len(_run_configs)
        dispersion_vals = []

        for step in range(S):
            x_sum = 0
            y_sum = 0
            for robot in range(len(_run_configs[step])):
                x_sum += _run_configs[step][robot][0]
                y_sum += _run_configs[step][robot][1]

            centroid = [x_sum / len(_run_configs[step]), y_sum / len(_run_configs[step])]

            dispersion_sum = 0
            for robot in range(len(_run_configs[step])):
                dispersion_sum += distArrs(_run_configs[step][robot], centroid)

            dispersion_vals.append(dispersion_sum)

        return dispersion_vals



    def convex_hull(self, _run_configs):
        """
        Takes as input an S x N x 3 run of position/orientation data and returns
        a 1 x S array of the value of the convex hull perimeter metric at each
        step.
        """
        S = len(_run_configs)
        convexhull_vals = []
        for step in range(S):
            robot_coords = []
            for robot in range(len(_run_configs[step])):
                robot_coords.append([ _run_configs[step][robot][0], _run_configs[step][robot][1] ])

            n = len(robot_coords)

            l = 0;
            for i in range(n):
              if (robot_coords[i][0] <= robot_coords[l][0]):
                l = i

            pointOnHull = robot_coords[l]
            firstPoint = pointOnHull
            maxThetaCandidates = []

            i = 0
            hull = []
            while True:
              hull.append(pointOnHull)
              if (i != 0):
                  robot_coords.remove(pointOnHull)

              endpoint = robot_coords[0]
              maxTheta = 0
              n = len(robot_coords)
              for j in range(n):
                if (i == 0):
                  beforeInitialPoint = [ hull[i][0], hull[i][1] - 10.0 ]
                  a = distArrs(hull[i], robot_coords[j]);
                  b = distArrs(beforeInitialPoint, robot_coords[j]);
                  c = distArrs(beforeInitialPoint, hull[i]);
                  if (a != 0 and b != 0 and c != 0):
                    if ( (c+a) == b ):
                      theta = math.pi
                    else:
                      theta = math.acos( ((a*a) + (c*c) - (b*b)) / (2 * a * c) )
                else:
                  a = distArrs(hull[i], robot_coords[j])
                  b = distArrs(hull[i-1], robot_coords[j])
                  c = distArrs(hull[i-1], hull[i])
                  if (a != 0 and b != 0 and c != 0):
                    if ( (c+a) == b ):
                      theta = math.pi
                    elif ( (a+b) == c ):
                      theta = 0.0
                    else:
                      if ( ( ((a*a) + (c*c) - (b*b)) / (2 * a * c) ) > 1 and
                           ( ((a*a) + (c*c) - (b*b)) / (2 * a * c) ) < 1.005 ):
                        theta = 0.0
                      elif ( ( ((a*a) + (c*c) - (b*b)) / (2 * a * c) ) < -1 and
                                ( ((a*a) + (c*c) - (b*b)) / (2 * a * c) ) > -1.005 ):
                        theta = math.pi
                      else:
                        theta = math.acos( ((a*a) + (c*c) - (b*b)) / (2 * a * c) )

                  if (i == 1):
                    if (robot_coords[j] == firstPoint):
                      theta = 0.0

                if (endpoint == pointOnHull):
                  endpoint = robot_coords[j]
                elif (theta > maxTheta):
                  maxTheta = theta
                  endpoint = robot_coords[j]
                  maxThetaCandidates.clear()
                  maxThetaCandidates.append(robot_coords[j])
                elif (theta == maxTheta):
                  maxThetaCandidates.append(robot_coords[j])
                  numCandidates = len(maxThetaCandidates)
                  maximumDist = 0.0
                  for z in range(numCandidates):
                    if ( distArrs(hull[i], maxThetaCandidates[z]) > maximumDist ):
                      endpoint = maxThetaCandidates[z]

              pointOnHull = endpoint
              i += 1
              if (endpoint == hull[0]):
                  break

            hn = len(hull)
            perimeter = 0.0
            for i in range(hn):
              if (i == hn-1):
                perimeter += distArrs(hull[i], hull[0])
              else:
                perimeter += distArrs(hull[i], hull[i+1])

            convexhull_vals.append(perimeter)

        return convexhull_vals




    def smallest_enclosing_disc_circum(self, _run_configs):
        """
        Takes as input an S x N x 3 run of position/orientation data and returns
        a 1 x S array of the value of the smallest enclosing disc circumference
        metric at each step.
        """
        S = len(_run_configs)
        sed_vals = []
        for step in range(S):
            robot_coords = []
            for robot in range(len(_run_configs[step])):
                robot_coords.append([ _run_configs[step][robot][0], _run_configs[step][robot][1] ])
            sed = welzl(robot_coords)
            sed_circum = sed.R * 2.0 * math.pi
            sed_vals.append(sed_circum)
        return sed_vals





    def cluster_fraction(self, _run_configs):
        """
        Takes as input an S x N x 3 run of position/orientation data and returns
        a 1 x S array of the value of the cluster fraction metric at each step.
        """
        S = len(_run_configs)
        cluster_fraction_vals = []
        for step in range(S):
            robot_coords = []
            for robot in range(len(_run_configs[step])):
                robot_coords.append([ _run_configs[step][robot][0], _run_configs[step][robot][1], 0 ])

            all_clusters = []
            for x in robot_coords:
                if(x[2] == 0):
                    current_cluster = []
                    DFS(x, robot_coords, current_cluster)
                    all_clusters.append(current_cluster)

            numInMaxCluster = len(all_clusters[0])
            for i in range(len(all_clusters)):
                if (len(all_clusters[i]) > numInMaxCluster):
                    numInMaxCluster = len(all_clusters[i])
            cluster_fraction_vals.append( numInMaxCluster / len(robot_coords) )
        return cluster_fraction_vals





    def animate(self, _run, _iter):
        """
        Animate the robots' movement over time on a 2D xy plane.
        """
        tqdm.write('Animating positions and movement of robots...')

        run_configs = self.runs_data[_run][_iter]
        S, N, _ = np.shape(run_configs)
        alpha = self.params[_run][2]
        time_step = self.params[_run][4]

        min = run_configs[0][0][0]
        max = run_configs[0][0][0]
        for step in range(S):
            for robot in range(N):
                if (run_configs[step][robot][0] < min):
                    min = run_configs[step][robot][0]
                if (run_configs[step][robot][0] > max):
                    max = run_configs[step][robot][0]
                if (run_configs[step][robot][1] < min):
                    min = run_configs[step][robot][1]
                if (run_configs[step][robot][1] > max):
                    max = run_configs[step][robot][1]

        fig, ax = plt.subplots(figsize=(5, 5), dpi=300)

        def roundup_hund(x):
            return int(math.ceil(x / 10.0)) * 10

        def rounddown_hund(x):
            return int(math.floor(x / 10.0)) * 10

        min_fig = rounddown_hund(min - 20)
        max_fig = roundup_hund(max + 20)

        ax.set(xlim=[min_fig, max_fig], ylim=[min_fig, max_fig])
        # plt.xticks(np.arange(min_fig, max_fig + 1, step=100))
        # plt.yticks(np.arange(min_fig, max_fig + 1, step=100))
        plt.title('Swarm Aggregation Simulation ('+ str(int(S/(1000/time_step))) + ' sec)')

        ims = []

        # cmap = np.vectorize(lambda x : cm.plasma(x))
        # colors = np.array(cmap(np.linspace(0.0, 1.0, num=N)))
        plasma = cm.get_cmap('plasma')
        colors = plasma(np.linspace(0,1,N))

        if (1000/time_step <= 60):
            frame_step = 1
            blit_frame_step = time_step
        else:
            frame_step = int(17/time_step)
            blit_frame_step = frame_step * time_step


        for step in range(0, S, frame_step):
            lines = []
            circles = []
            sight_cones = []
            for robot in range(len(run_configs[step])):
                x1 = run_configs[step][robot][0]
                y1 = run_configs[step][robot][1]
                x2 = run_configs[step][robot][0] + 100000*math.cos(run_configs[step][robot][2] + (math.pi/2))
                y2 = run_configs[step][robot][1] + 100000*math.sin(run_configs[step][robot][2] + (math.pi/2))
                lines.append([(x1,y1), (x2,y2)])
                circles.append(plt.Circle( (run_configs[step][robot][0], run_configs[step][robot][1]), radius=3.7, linewidth=0, color=colors[robot] ) )
                cone_x1=100000
                cone_x2=100000
                if ( ( run_configs[step][robot][2]+(math.pi/2) +(alpha/2) ) %(2*math.pi) > (math.pi/2) and ( run_configs[step][robot][2]+(math.pi/2) +(alpha/2) ) %(2*math.pi) < (3*math.pi/2)):
                    cone_x1 = cone_x1 * (-1)
                if ( ( run_configs[step][robot][2]+(math.pi/2) -(alpha/2) ) %(2*math.pi) > (math.pi/2) and ( run_configs[step][robot][2]+(math.pi/2) -(alpha/2) ) %(2*math.pi) < (3*math.pi/2)):
                    cone_x2 = cone_x2 * (-1)
                cone_y1 = math.tan(run_configs[step][robot][2]+(math.pi/2) +(alpha/2))*(cone_x1-x1) + y1
                cone_y2 = math.tan(run_configs[step][robot][2]+(math.pi/2) -(alpha/2))*(cone_x2-x1) + y1
                if ( (run_configs[step][robot][2]+(math.pi/2) +(alpha/2)) %(2*math.pi) == math.pi/2):
                    cone_y1 = y1 + 100000
                    cone_x1 = x1
                elif ( (run_configs[step][robot][2]+(math.pi/2) +(alpha/2)) %(2*math.pi) == 3*math.pi/2):
                    cone_y1 = y1 - 100000
                    cone_x1 = x1
                if ( (run_configs[step][robot][2]+(math.pi/2) -(alpha/2)) %(2*math.pi) == math.pi/2):
                    cone_y2 = y1 + 100000
                    cone_x2 = x1
                elif ( (run_configs[step][robot][2]+(math.pi/2) -(alpha/2)) %(2*math.pi) == 3*math.pi/2):
                    cone_y2 = y1 - 100000
                    cone_x2 = x1
                sight_cones.append(plt.Polygon([[x1,y1], [cone_x1,cone_y1], [cone_x2,cone_y2]], color=colors[robot], alpha=0.3))
            # col = LineCollection(lines, colors=colors.transpose(), linewidths=.5)
            c = PatchCollection(circles, match_original=True)
            col = LineCollection(lines, colors=colors, linewidths=.5, alpha=.75)
            cones = PatchCollection(sight_cones, match_original=True)
            # ims.append([ ax.add_collection(col), ax.scatter([data[0] for data in run_configs[step]], [data[1] for data in run_configs[step]], c='C0', s=54.76) ])
            ims.append([ ax.add_collection(cones), ax.add_collection(col), ax.add_collection(c) ])

        ani = matplotlib.animation.ArtistAnimation(fig, ims, interval=blit_frame_step, blit=True)
        ani.save('animations/' + self.fname + '_run' + str(_run) + '_iter' + str(_iter) + 'animation.mp4')




    # EXP 0: time evolutions; plot metrics over time for a single iter of a single run.
    def plot_exp0(self, _run, _iter):
        tqdm.write('Plotting metrics over time...')

        run_configs = self.runs_data[_run][_iter]

        time_step = self.params[_run][4]

        # _metrics = ['dispersion', 'convex_hull', 'smallest_enclosing_disc_circum', 'cluster_fraction']
        # _metrics = ['dispersion', 'convex_hull', 'cluster_fraction']
        _metrics = ['dispersion']

        for i in range(len(_metrics)):
            fig, ax = plt.subplots()
            if _metrics[i] == 'dispersion':
                y = self.dispersion(run_configs)
                metric_str_name = 'Dispersion'
            elif _metrics[i] == 'convex_hull':
                y = self.convex_hull(run_configs)
                metric_str_name = 'Convex Hull Perim.'
            elif _metrics[i] == 'smallest_enclosing_disc_circum':
                y = self.smallest_enclosing_disc_circum(run_configs)
                metric_str_name = 'Smallest Enclosing Disc Circumf.'
            else:   # _metrics[i] == 'cluster_fraction'
                y = self.cluster_fraction(run_configs)
                metric_str_name = 'Cluster Fraction'
            ax.plot(np.arange(0, time_step*len(y)/1000, time_step/1000), y, color='C0')
            ax.set(xlabel='Time (sec)', ylabel=metric_str_name)
            ax.grid()
            plt.tight_layout()
            fig.savefig('figs/' + self.fname + '_' + _metrics[i] + '_run' + str(_run) + '_iter' + str(_iter) + '.png', dpi=300)
            plt.close()





    def plot_expsymm(self, _run, _iter):
        tqdm.write('Plotting metrics over time...')

        run_configs = self.runs_data[_run][_iter]
        fig, ax = plt.subplots()
        y = self.dispersion(run_configs)


        abs_min=y[0]
        for i in range(len(y)):
            if(y[i]<abs_min):
                abs_min=y[i]
        print(abs_min)


        ax.plot(np.arange(0, len(y)/2000, 1/2000), y)
        ax.set(xlabel='Time (sec)', ylabel='Dispersion')
        ax.grid()
        plt.tight_layout()
        fig.savefig('figs/' + self.fname + '_symm_dispersion.png', dpi=300)
        plt.close()





    # Error Prob
    def plot_exp1_errorprob(self):
        print('Plotting metrics over time...')

        from operator import add
        from operator import sub

        def aggregate_stats(data):
            """
            Calculates the average and standard deviation of data by row.
            Inputs:  A list of lists where each row represents data for a different value of the independent variable and each column represents values obtained from repeated experiment runs.
            Returns: A dictionary of two lists, one representing each row's average and the other representing its standard deviation.
            """
            return {'ave' : [stats.mean(row) for row in data], \
                    'stddev' : [stats.pstdev(row) for row in data]}

        def plot_errortube(ax, x, y, yerr, color):
            """
            Plots a translucent error tube centered at y with error sizes yerr with the
            specified color.
            Inputs:  pyplot axes to draw on, lists of x and y coordinates, a list of
                     error magnitudes at each x coordinate, and a hex color to draw the
                     error tube in.
            Returns: None.
            """
            errs_lower = list(map(sub, y, yerr))
            errs_upper = list(map(add, y, yerr))
            ax.fill_between(x, errs_lower, errs_upper, alpha=-0.5, facecolor=color)


        params = np.copy(self.params)
        ten_data = []
        twenfi_data = []
        fifty_data = []
        onehund_data = []

        for run in range(len(params)):
            if (params[run][0] == 10) :
                new_row = []
                for iter in range(self.iters):
                    new_row.append(self.runs_data[run][iter]/2000)
                ten_data.append(new_row)
            if (params[run][0] == 25) :
                new_row = []
                for iter in range(self.iters):
                    new_row.append(self.runs_data[run][iter]/2000)
                twenfi_data.append(new_row)
            elif (params[run][0] == 50) :
                new_row = []
                for iter in range(self.iters):
                    new_row.append(self.runs_data[run][iter]/2000)
                fifty_data.append(new_row)
            elif (params[run][0] == 100) :
                new_row = []
                for iter in range(self.iters):
                    new_row.append(self.runs_data[run][iter]/2000)
                onehund_data.append(new_row)

        scatter_twenfi_data = []
        scatter_x_twenfi = []
        for run in range(len(twenfi_data)):
            x = 0.0 + run*.0025
            for iter in range(len(twenfi_data[run])):
                scatter_twenfi_data.append(twenfi_data[run][iter])
                scatter_x_twenfi.append(x)

        x = [0.0 + (i*.0025) for i in range(len(twenfi_data))]

        ten_data = aggregate_stats(ten_data)
        twenfi_data = aggregate_stats(twenfi_data)
        fifty_data = aggregate_stats(fifty_data)
        onehund_data = aggregate_stats(onehund_data)

        # cmap = np.vectorize(lambda x : cm.plasma(x))
        # colors = np.array(cmap(np.linspace(0.0, 1.0, num=4)))
        plasma = cm.get_cmap('plasma')
        colors = plasma(np.linspace(0,1,4))

        fig,ax = plt.subplots()

        ax.plot(x, np.full(len(x), 300), color='#000000', linestyle='dashed', alpha=0.6, zorder=2)

        # ax.plot(x, ten_data['ave'], color=colors[0], label='10 robots')
        # plot_errortube(ax, x, ten_data['ave'], ten_data['stddev'], to_hex(colors[0]))

        ax.plot(x, twenfi_data['ave'], color=colors[0], label='25 robots', zorder=4, linewidth=2)
        ax.scatter(scatter_x_twenfi, scatter_twenfi_data, color=colors[0], alpha=0.5, zorder=3)

        # ax.plot(x, fifty_data['ave'], color=colors[2], label='50 robots')
        # plot_errortube(ax, x, fifty_data['ave'], fifty_data['stddev'], to_hex(colors[2]))
        #
        # ax.plot(x, onehund_data['ave'], color=colors[3], label='100 robots')
        # plot_errortube(ax, x, onehund_data['ave'], onehund_data['stddev'], to_hex(colors[3]))

        ax.set(title='Effect of Error Probability Noise on Runtime', xlabel='Error Probability', ylabel='Runtime (Seconds)')
        ax.legend(loc='lower right')
        ax.grid()


        fig.savefig('figs/' + self.fname + '_scatter.png', dpi=300)
        plt.close()





    # Motion
    def plot_exp1_motion(self):
        print('Plotting metrics over time...')

        from operator import add
        from operator import sub

        def aggregate_stats(data):
            """
            Calculates the average and standard deviation of data by row.
            Inputs:  A list of lists where each row represents data for a different value of the independent variable and each column represents values obtained from repeated experiment runs.
            Returns: A dictionary of two lists, one representing each row's average and the other representing its standard deviation.
            """
            return {'ave' : [stats.mean(row) for row in data], \
                    'stddev' : [stats.pstdev(row) for row in data]}

        def plot_errortube(ax, x, y, yerr, color):
            """
            Plots a translucent error tube centered at y with error sizes yerr with the
            specified color.
            Inputs:  pyplot axes to draw on, lists of x and y coordinates, a list of
                     error magnitudes at each x coordinate, and a hex color to draw the
                     error tube in.
            Returns: None.
            """
            errs_lower = list(map(sub, y, yerr))
            errs_upper = list(map(add, y, yerr))
            ax.fill_between(x, errs_lower, errs_upper, alpha=-0.5, facecolor=color)


        params = np.copy(self.params)
        ten_data = []
        twenfi_data = []
        fifty_data = []
        onehund_data = []

        for run in range(len(params)):
            if (params[run][0] == 10) :
                new_row = []
                for iter in range(self.iters):
                    new_row.append(self.runs_data[run][iter]/2000)
                ten_data.append(new_row)
            if (params[run][0] == 25) :
                new_row = []
                for iter in range(self.iters):
                    new_row.append(self.runs_data[run][iter]/2000)
                twenfi_data.append(new_row)
            elif (params[run][0] == 50) :
                new_row = []
                for iter in range(self.iters):
                    new_row.append(self.runs_data[run][iter]/2000)
                fifty_data.append(new_row)
            elif (params[run][0] == 100) :
                new_row = []
                for iter in range(self.iters):
                    new_row.append(self.runs_data[run][iter]/2000)
                onehund_data.append(new_row)


        scatter_ten_data = []
        scatter_x_ten = []
        for run in range(len(ten_data)):
            x = 0.0 + run*.25
            for iter in range(len(ten_data[run])):
                scatter_ten_data.append(ten_data[run][iter])
                scatter_x_ten.append(x)

        x = [0.0 + (i*.25) for i in range(len(ten_data))]

        ten_data = aggregate_stats(ten_data)
        twenfi_data = aggregate_stats(twenfi_data)
        fifty_data = aggregate_stats(fifty_data)
        onehund_data = aggregate_stats(onehund_data)

        plasma = cm.get_cmap('plasma')
        colors = plasma(np.linspace(0,1,4))

        fig,ax = plt.subplots()


        ax.plot(x, np.full(len(x), 300), color='#000000', linestyle='dashed', alpha=0.6, zorder=2)

        ax.plot(x, ten_data['ave'], color=colors[0], label='10 robots', zorder=4, linewidth=2)
        ax.scatter(scatter_x_ten, scatter_ten_data, color=colors[0], alpha=0.5, zorder=3)

        # ax.plot(x, twenfi_data['ave'], color=colors[1], label='25 robots')
        # plot_errortube(ax, x, twenfi_data['ave'], twenfi_data['stddev'], to_hex(colors[1]))
        #
        # ax.plot(x, fifty_data['ave'], color=colors[2], label='50 robots')
        # plot_errortube(ax, x, fifty_data['ave'], fifty_data['stddev'], to_hex(colors[2]))
        #
        # ax.plot(x, onehund_data['ave'], color=colors[3], label='100 robots')
        # plot_errortube(ax, x, onehund_data['ave'], onehund_data['stddev'], to_hex(colors[3]))

        ax.set(title='Effect of Motion Noise on Runtime', xlabel='Maximum Motion Noise Force (N)', ylabel='Runtime (Seconds)')
        ax.legend(loc='lower right')
        ax.grid()


        fig.savefig('figs/' + self.fname + '_scatter.png', dpi=300)
        plt.close()






    def plot_expsymm2(self):
        tqdm.write('Plotting metrics over time...')

        three_robots = self.dispersion(self.runs_data[0][0])
        x3=len(three_robots)
        five_robots = self.dispersion(self.runs_data[1][0])
        ten_robots = self.dispersion(self.runs_data[2][0])

        plasma = cm.get_cmap('plasma')
        colors = plasma(np.linspace(0,1,4))

        fig,ax = plt.subplots()

        ax.plot(np.arange(0, len(three_robots)/2000, 1/2000), three_robots, color=colors[0], label='3 robots')
        ax.plot(np.arange(0, len(three_robots)/2000, 1/2000), np.full(len(three_robots), 12.817176), color=colors[0], linestyle='dashed')

        ax.plot(np.arange(0, len(five_robots)/2000, 1/2000), five_robots, color=colors[1], label='5 robots')
        ax.plot(np.arange(0, len(five_robots)/2000, 1/2000), np.full(len(five_robots), 28.8987147321), color=colors[1], linestyle='dashed')

        ax.plot(np.arange(0, len(ten_robots)/2000, 1/2000), ten_robots, color=colors[2], label='10 robots')
        ax.plot(np.arange(0, len(ten_robots)/2000, 1/2000), np.full(len(ten_robots), 81.57429538), color=colors[2], linestyle='dashed')

        ax.set(title='Symmetric Initial Configurations', xlabel='Time (sec)', ylabel='Dispersion')
        ax.legend(loc='upper right')
        ax.grid()
        plt.tight_layout()
        fig.savefig('figs/' + self.fname + '_symm_dispersion.png', dpi=300)
        plt.close()








# Load data - noise
def plot_load_noise(runtimes, filename):
    print('Plotting metrics over time...')

    from operator import add
    from operator import sub

    def aggregate_stats(data):
        """
        Calculates the average and standard deviation of data by row.
        Inputs:  A list of lists where each row represents data for a different value of the independent variable and each column represents values obtained from repeated experiment runs.
        Returns: A dictionary of two lists, one representing each row's average and the other representing its standard deviation.
        """
        return {'ave' : [stats.mean(row) for row in data], \
                'stddev' : [stats.pstdev(row) for row in data]}

    def plot_errortube(ax, x, y, yerr, color):
        """
        Plots a translucent error tube centered at y with error sizes yerr with the
        specified color.
        Inputs:  pyplot axes to draw on, lists of x and y coordinates, a list of
                 error magnitudes at each x coordinate, and a hex color to draw the
                 error tube in.
        Returns: None.
        """
        errs_lower = list(map(sub, y, yerr))
        errs_upper = list(map(add, y, yerr))
        ax.fill_between(x, errs_lower, errs_upper, alpha=-0.5, facecolor=color)


    ten_data = []
    for run in range(len(runtimes[0])):
        new_row = []
        for iter in range(len(runtimes[0][run])):
            new_row.append(runtimes[0][run][iter]/2000)
        ten_data.append(new_row)

    twenfi_data = []
    for run in range(len(runtimes[1])):
        new_row = []
        for iter in range(len(runtimes[1][run])):
            new_row.append(runtimes[1][run][iter]/2000)
        twenfi_data.append(new_row)

    x = [0.0 + (i*(math.pi/24)) for i in range(len(ten_data))]

    ten_data = aggregate_stats(ten_data)
    twenfi_data = aggregate_stats(twenfi_data)

    plasma = cm.get_cmap('plasma')
    colors = plasma(np.linspace(0,1,4))

    fig,ax = plt.subplots()

    ax.plot(x, ten_data['ave'], color=colors[0], label='10 robots')
    plot_errortube(ax, x, ten_data['ave'], ten_data['stddev'], to_hex(colors[0]))

    ax.plot(x, twenfi_data['ave'], color=colors[1], label='25 robots')
    plot_errortube(ax, x, twenfi_data['ave'], twenfi_data['stddev'], to_hex(colors[1]))


    ax.set(title='Effect of Sight Sensor Size on Runtime', xlabel='Cone-of-sight Sensor Size (rad)', ylabel='Runtime (Seconds)')
    ax.legend(loc='upper right')
    ax.grid()


    fig.savefig('figs/' + filename + '_sightcone.png', dpi=300)
    plt.close()







# Load data - cone-of-sight sensor size
def plot_sightcone(runtimes, filename):
    print('Plotting metrics over time...')

    from operator import add
    from operator import sub

    def aggregate_stats(data):
        """
        Calculates the average and standard deviation of data by row.
        Inputs:  A list of lists where each row represents data for a different value of the independent variable and each column represents values obtained from repeated experiment runs.
        Returns: A dictionary of two lists, one representing each row's average and the other representing its standard deviation.
        """
        return {'ave' : [stats.mean(row) for row in data], \
                'stddev' : [stats.pstdev(row) for row in data]}

    def plot_errortube(ax, x, y, yerr, color):
        """
        Plots a translucent error tube centered at y with error sizes yerr with the
        specified color.
        Inputs:  pyplot axes to draw on, lists of x and y coordinates, a list of
                 error magnitudes at each x coordinate, and a hex color to draw the
                 error tube in.
        Returns: None.
        """
        errs_lower = list(map(sub, y, yerr))
        errs_upper = list(map(add, y, yerr))
        ax.fill_between(x, errs_lower, errs_upper, alpha=-0.5, facecolor=color)


    ten_data = []
    x_ten = []
    for run in range(len(runtimes[0])):
        x = 0.0 + run*(math.pi/24)
        for iter in range(len(runtimes[0][run])):
            ten_data.append(runtimes[0][run][iter]/2000)
            x_ten.append(x)

    twenfi_data = []
    x_twenfi = []
    for run in range(len(runtimes[1])):
        x = 0.0 + run*(math.pi/24)
        for iter in range(len(runtimes[1][run])):
            twenfi_data.append(runtimes[1][run][iter]/2000)
            x_twenfi.append(x)
    for run in range(len(runtimes[2])):
        x = 0.0 + run*(math.pi/24)
        for iter in range(len(runtimes[2][run])):
            twenfi_data.append(runtimes[2][run][iter]/2000)
            x_twenfi.append(x)

    main_x = []
    ten_avg = []
    twenfi_avg = []
    for i in range(25):
        main_x.append(i * (math.pi/24))
        ten_avg.append([])
        twenfi_avg.append([])

    for i in range(len(ten_data)):
        x = int( (x_ten[i] / math.pi * 24) + .001)
        ten_avg[x].append(ten_data[i])

    for i in range(len(twenfi_data)):
        x = int( (x_twenfi[i] / math.pi * 24) + .001)
        twenfi_avg[x].append(twenfi_data[i])

    ten_avg = aggregate_stats(ten_avg)
    twenfi_avg = aggregate_stats(twenfi_avg)


    plasma = cm.get_cmap('plasma')
    colors = plasma(np.linspace(0,1,4))

    fig,ax = plt.subplots()

    ax.plot(x_ten, np.full(len(x_ten), 300), color='#000000', linestyle='dashed', alpha=0.6, zorder=2)

    ax.scatter(x_ten, ten_data, color=colors[0], alpha=0.4, label='10 robots', zorder=3)
    ax.plot(main_x, ten_avg['ave'], color=colors[0], zorder=5)

    ax.scatter(x_twenfi, twenfi_data, color=colors[2], alpha=0.4, label='25 robots', zorder=4)
    ax.plot(main_x, twenfi_avg['ave'], color=colors[2], zorder=6)


    ax.set(title='Effect of Sight Sensor Size on Runtime', xlabel='Cone-of-sight Sensor Size (rad)', ylabel='Runtime (Seconds)')
    ax.legend(loc='upper right')
    ax.grid()

    plt.ylim((-5, 315))

    fig.savefig('figs/' + filename + '_sightcone.png', dpi=300)
    plt.close()







def load(filename):
    """
    Loads the experiment data from a previously saved pickle file.
    """
    tqdm.write('Loading Experiment ' + filename + '...')
    file = open('data/'+ filename, 'rb')
    data = pickle.load(file)
    file.close()
    return data









# Evidence for symmetric counterexample
def expsymm(_seed=None):
    params = {'N' : [3, 5, 10], 'S' : [4 * 60 * 2*1000], 'T' : [0], 'noise_p' : [0]}
    exp = Experiment(_id='symm', _params=params, _init='symmetric', _noise='errorprob', _stopping=False, _iters=1, _seed=_seed)
    exp.run()
    exp.save()
    exp.plot_expsymm2()
    print(exp.fname)

# Time evolutions
def exp0(_seed=None):
    params = {'N' : [100], 'S' : [300 * 2*1000], 'T' : [0], 'noise_p' : [0.00]}
    exp = Experiment(_id='0', _params=params, _init='random', _noise='errorprob', _stopping=False, _savehistory=True, _iters=1, _seed=_seed)
    exp.run()
    exp.save()
    print(exp.fname)
    exp.plot_exp0(0, 0)
    exp.animate(0, 0)

def exp_test(_seed=None):
    params = {'N' : [10, 25, 50], 'S' : [0], 'T' : [0], 'noise_p' : [1], 'time_step' : [10]}
    exp = Experiment(_id='test', _params=params, _init='random', _noise='motion', _stopping=True, _savehistory=False, _iters=2, _seed=_seed)
    exp.run()
    exp.save()
    print(exp.fname)
    # print(exp.runs_data)
    # exp.plot_exp0(0,0)
    # exp = load('exptest_458736_172917691849.pkl')
    # for run in range(len(exp.params)):
    #     exp.animate(run,0)

# Error probability noise
def exp1_errorprob(_seed=None):
    params = {'N' : [10, 25, 50, 100], 'S' : [0], 'T' : [0], 'noise_p' : np.arange(0, .1000001, .002)}
    exp = Experiment(_id='1errorprob', _params=params, _init='random', _noise='errorprob', _stopping=True, _savehistory=False, _iters=20, _seed=_seed)
    exp.run()
    exp.save()
    print(exp.fname)
    # exp.plot_exp1_errorprob()

# 'Motion' noise
def exp1_motion(_seed=None):
    # params = {'N' : [10, 25, 50, 100], 'S' : [0], 'T' : [0], 'noise_p' : np.arange(0, 10.0001, .2)}
    params = {'N' : [10, 25], 'S' : [0], 'T' : [0], 'noise_p' : np.arange(0, 10.0001, .4)}
    exp = Experiment(_id='1motion', _params=params, _init='random', _noise='motion', _stopping=True, _savehistory=False, _iters=10, _seed=_seed)
    exp.run()
    exp.save()
    print(exp.fname)
    # exp.plot_exp1_motion()

    # exp.plot_exp0(0,0)
    # run=0
    # while(run <= 2):
    #     exp.plot_exp0(run,0)
    #     exp.plot_exp0(run,1)
    #     run+=1

# Cone of sight vs. line of sight
def exp2(_seed=None):
    params = {'N' : [50], 'S' : [0], 'T' : np.arange(0, math.pi+.00001, math.pi/24), 'noise_p' : [0.0]}
    exp = Experiment(_id='2', _params=params, _init='random', _noise='errorprob', _stopping=True, _savehistory=False, _iters=10, _seed=_seed)
    exp.run()
    exp.save()
    print(exp.fname)
    # exp.plot_exp2()

def expload(_seed=None):
    # # filename='exp1errorprob_9123735_111435885763.pkl'
    # data = [load('exp2_8124652312_115657246372.pkl').runs_data, load('exp2_947562543_182100626627.pkl').runs_data, load('exp2_4239751634_235948089486.pkl').runs_data]
    # # print(data.runs_data)
    # # print(data.params[1])
    # plot_sightcone(data, 'exp2_8124652312_115657246372')

    exp = load('exp1errorprob_9123735_111435885763.pkl')
    print(exp.params[1])
    exp.plot_exp1_errorprob()

def expdeadlock(_seed=None):
    params = {'N' : [17], 'S' : [60 * 2*1000], 'T' : [0], 'noise_p' : [0.0]}
    exp = Experiment(_id='deadlock', _params=params, _init='custom', _noise='errorprob', _stopping=False, _savehistory=True, _iters=1, _seed=_seed)
    exp.run()
    exp.save()
    print(exp.fname)
    exp.plot_exp0(0,0)
    exp.animate(0,0)

def exptimestep(_seed=None):
    params = {'N' : [25], 'S' : [60 * 3], 'T' : [0], 'noise_p' : [0.0], 'time_step' : [1, 2, 5, 10, 50, 100]}
    exp = Experiment(_id='timestep', _params=params, _init='random', _noise='errorprob', _stopping=False, _savehistory=True, _iters=1, _seed=_seed)
    exp.run()
    exp.save()
    print(exp.fname)
    for i in range(6):
        exp.plot_exp0(i,0)
        exp.animate(i,0)




if __name__ == '__main__':
    # Parse command line arguments.
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-E', '--exps', type=str, nargs='+', required=True, \
                        help='Indices of experiments to run')
    parser.add_argument('-R', '--rand_seed', type=int, default=None, \
                        help='Seed for random number generation')
    args = parser.parse_args()


    # Run selected experiments.
    exps = {'symm' : expsymm, '0' : exp0, 'test' : exp_test, '1errorprob' : exp1_errorprob, '1motion' : exp1_motion, '2' : exp2, 'load' : expload, 'deadlock' : expdeadlock, 'timestep' : exptimestep}
    for id in args.exps:
        exps[id](args.rand_seed)
