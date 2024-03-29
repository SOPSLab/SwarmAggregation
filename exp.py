# Project:     SwarmAggregation
# Filename:    exp.py
# Authors:     Joshua J. Daymude (jdaymude@asu.edu) and Noble C. Harasha
#              (nharasha@mit.edu).

"""
exp: A flexible, unifying framework for defining and running experiments for
     swarm aggregation.
"""

import argparse
from aggregation import aggregation, ideal
from itertools import product
from math import sin, cos, hypot, ceil
from matplotlib.animation import FFMpegWriter, ArtistAnimation
import matplotlib.cm as cm
from matplotlib.collections import LineCollection, PatchCollection, PolyCollection
import matplotlib.pyplot as plt
from metrics import *
import numpy as np
import pickle
from tqdm import tqdm


class Experiment(object):
    """
    A flexible, unifying framework for experiments.
    """

    def __init__(self, id, params={}, iters=1, savehist=True, seed=None):
        """
        Inputs:
        - id (str): identifier for the experiment
        - params (dict): the full parameter set for the simulation runs
            {
              'N' : [int > 0] number of robots,
              'R' : [float > 0] radius of rotation (m),
              'r' : [float > 0] radius of a robot (m),
              'm' : [float > 0] mass of a robot (kg),
              'w0' : [float] rot. speed of a robot about its center (rad/s),
              'w1' : [float] rot. speed of a robot in place (rad/s),
              'sensor' : [0 <= float <= pi] size of the sight sensor (rad),
              'noise' : [(str, float)] either ('err', p) for error probability
                        with probability p or ('mot', f) for motion noise with
                        maximum force f (N),
              'time' : [float > 0] wall-clock duration of simulation (s),
              'step' : [float > 0] wall-clock duration of a time step (s),
              'stop' : [float >= 0] if not None, simulation stops if system's
                       dispersion is within stop% of the ideal value,
              'init' : ['rand', 'symm'] initialization mode
            }
        - iters (int): the number of iterated runs for each parameter setting
        - savehist (bool): True if a run's history should be saved
        - seed (int): random seed
        """
        # Unpack singular parameters.
        self.id, self.iters, self.savehist, self.seed, = id, iters, savehist, seed

        # Unpack aggregation parameters.
        defaults = {'N' : [100], 'R' : [0.1445], 'r' : [0.037], 'm' : [0.125], \
                    'w0' : [-0.75], 'w1' : [-5.02], 'sensor' : [0], \
                    'noise' : [('err', 0)], 'time' : [300], 'step' : [0.005], \
                    'stop' : [None], 'init' : ['rand']}
        plist = [params[p] if p in params else defaults[p] for p in defaults]
        self.params = list(product(*plist))

        # Set up data and results filenames.
        self.fname = 'exp_{}_{}'.format(self.id, self.seed)

        # Instantiate a list to hold runs data. This data will have shape
        # A x B x [S x N x 3, 1] where A is the number of runs (i.e., unique
        # parameter combinations), B is the number of iterations per run, S is
        # the number of time steps simulated, N is the number of robots, and 3
        # represents each robot's X/Y/Theta data.
        self.runs_data = [[] for p in self.params]


    def run(self):
        """
        Run this experiment according to the input parameters.
        """
        tqdm.write('Running Experiment ' + self.id + '...')

        # Set up random seeds for iterated runs.
        rng = np.random.default_rng(self.seed)
        run_seeds = rng.integers(0, 2**32, size=self.iters)

        # For each parameter combination, do iterated runs of aggregation.
        silent = len(self.params) > 1 or self.iters > 1
        for i, param in enumerate(tqdm(self.params, desc='Simulating runs')):
            N, R, r, m, w0, w1, sensor, noise, time, step, stop, init = param
            for seed in tqdm(run_seeds, desc='Iterating run', \
                             leave=bool(i == len(self.params) - 1)):
                run_data = aggregation(N, R, r, m, w0, w1, sensor, noise, time,\
                                       step, stop, init, seed, silent)
                if not self.savehist:
                    # Only save the final configuration.
                    history, final = run_data
                    self.runs_data[i].append((np.copy(history[final-1]), final))
                else:
                    # Save the entire configuration history.
                    self.runs_data[i].append(run_data)


    def save(self):
        """
        Saves this experiment, including all parameters and run data, to a file
        named according to the experiment's ID and seed.
        """
        tqdm.write('Saving Experiment ' + self.id + '...')
        with open('data/' + self.fname + '.pkl', 'wb') as f:
            pickle.dump(self, f)


    def plot_evo(self, runs, iters, metrics=['sed', 'hull', 'disp', 'clus'], \
                 labels=None, title='', anno=''):
        """
        Takes indices of either (i) one run and multiple iterations or (ii) one
        iteration of multiple runs and plots the given metrics against time.
        """
        tqdm.write('Plotting metrics over time...')

        # Sanity checks and setup. Assumes N, r, time, and step are static.
        assert self.savehist, 'ERROR: No history to calculate metrics per step'
        assert len(runs) == 1 or len(iters) == 1, 'ERROR: One run or one iter'
        runits = [i for i in product(runs, iters)]

        # Set up colors.
        cmap = np.vectorize(lambda x : cm.inferno(x))
        c = np.array(cmap(np.linspace(0, 1, len(runits) + 2))).T

        # Plot metrics over time for each run/iteration.
        names = {'sed'  : 'Smallest Enclosing Disc Circumference', \
                 'hull' : 'Convex Hull Perimeter', \
                 'disp' : 'Dispersion', \
                 'clus' : 'Cluster Fraction'}
        for metric in metrics:
            fig, ax = plt.subplots()
            for i, runit in enumerate(tqdm(runits)):
                # Plot the given metric over time.
                N, r, time, step = [self.params[runit[0]][j] for j in [0,2,8,9]]
                configs, final = self.runs_data[runit[0]][runit[1]]
                x = np.arange(0, time + step, step)[:final]
                y = []
                for config in tqdm(configs, desc='Calculating '+names[metric]):
                    if metric == 'sed':
                        y.append(sed_circumference(config))
                    elif metric == 'hull':
                        y.append(hull_perimeter(config))
                    elif metric == 'disp':
                        y.append(dispersion(config))
                    else:  # metric == 'clus'
                        y.append(cluster_fraction(config, r))
                if labels != None:
                    ax.plot(x, y, color=c[i+1], label=labels[i], zorder=4)
                else:
                    ax.plot(x, y, color=c[i+1], zorder=4)

                # Plot the minimum value for this metric as a dashed line.
                if metric == 'sed':
                    metric_min = sed_circumference(ideal(N, r))
                elif metric == 'hull':
                    metric_min = hull_perimeter(ideal(N, r))
                elif metric == 'disp':
                    metric_min = dispersion(ideal(N, r))
                else:  # metric == 'clus'
                    metric_min = cluster_fraction(ideal(N, r), r)
                ax.plot(x, np.full(len(x), metric_min), color=c[i+1], \
                        linestyle='dashed', zorder=3)

            # Save figure.
            ax.set(title=title, xlabel='Time (s)', ylabel=names[metric])
            ax.set_ylim(bottom=0)
            ax.grid()
            if labels != None:
                ax.legend(loc='upper right')
            plt.tight_layout()
            fig.savefig('figs/' + self.fname + '_' + metric + anno + '.png', \
                        dpi=300)
            plt.close()


    def plot_aggtime(self, N, ps, plabel, title='', anno=''):
        """
        Plots final and average time to aggregation per parameter value per
        number of robots. Assumes that the only parameters that are varied are
        the number of robots (N) and one non-time related parameter.
        """
        tqdm.write('Plotting average time to aggregation...')

        # Set up figure and colors.
        fig, ax = plt.subplots()
        cmap = np.vectorize(lambda x : cm.inferno(x))
        c = np.array(cmap(np.linspace(0, 1, len(N) + 2))).T

        # Plot simulation time cutoff as a dashed line.
        time, step = self.params[0][8], self.params[0][9]
        ax.plot(ps, np.full(len(ps), time), color='k', linestyle='dashed')

        # Plot iteration times as a scatter plot and averages as lines.
        for i, ni in enumerate(N):
            xs, ys, aves = [], [], []
            for j, run in enumerate(self.runs_data[i*len(ps):(i+1)*len(ps)]):
                agg_times = []
                for iter in run:
                    xs.append(ps[j])
                    agg_times.append(iter[1] * step)
                ys += agg_times
                aves.append(np.mean(agg_times))
            ax.scatter(xs, ys, color=c[i+1], s=15, alpha=0.4)
            ax.plot(ps, aves, color=c[i+1], label='{} robots'.format(ni))

        # Save figure.
        ax.set(title=title, xlabel=plabel, ylabel='Aggregation Time (s)')
        ax.set_ylim(bottom=0)
        ax.grid()
        ax.legend(loc='upper left')
        plt.tight_layout()
        fig.savefig('figs/' + self.fname + '_aggtime' + anno + '.png', dpi=300)
        plt.close()


    def animate(self, run, iter, frame=25, anno=''):
        """
        Animate the robots' movement over time.
        """
        tqdm.write('Animating robots\' movement...')

        # Check that a configuration history exists.
        assert self.savehist, 'ERROR: No history to animate'

        # Check that the desired frame rate is valid.
        assert frame > 0, 'ERROR: Frame rate must be positive value'

        # Get data and parameters.
        configs, final = self.runs_data[run][iter]
        N, r, sensor, time, step = [self.params[run][i] for i in [0,2,6,8,9]]

        # Set up plot.
        fig, ax = plt.subplots(figsize=(5,5), dpi=300)
        all_xy = configs[:,:,:2].flatten()
        fig_min, fig_max = np.min(all_xy) - r, np.max(all_xy) + r
        ax.set(xlim=[fig_min, fig_max], ylim=[fig_min, fig_max])

        # Set up colors for the various robots.
        cmap = np.vectorize(lambda x : cm.inferno(x))
        c = np.array(cmap(np.linspace(0, 0.9, N))).T

        # Set up frame rate to target at most 'frame' fps in real time.
        frame_step = 1 if step >= 1 / frame else ceil(1 / frame / step)
        interval = (step * frame_step) * 1000  # ms

        ims = []
        max_dist = hypot(*np.full(2, fig_max-fig_min))
        for s in tqdm(np.arange(0, min(len(configs), final), frame_step)):
            title = plt.text(1.0, 1.02, '{:.2f}s of {}s'.format(s*step, time), \
                             ha='right', va='bottom', transform=ax.transAxes)
            robots, lines, cones = [], [], []
            for i in range(N):
                xy, theta = configs[s][i][:2], configs[s][i][2]
                sensor_xy = xy + np.array([r * cos(theta), r * sin(theta)])

                # Add this robot's circle artist.
                robots.append(plt.Circle(xy, radius=r, linewidth=0, color=c[i]))

                # Add this robot's sight sensor direction artist.
                vec = max_dist * np.array([cos(theta), sin(theta)])
                lines.append([sensor_xy, sensor_xy + vec])

                # Add this robot's cone-of-sight polygon artist.
                if sensor > 0:
                    cw, ccw = theta - sensor / 2, theta + sensor / 2
                    vec_cw = max_dist * np.array([cos(cw), sin(cw)])
                    vec_ccw = max_dist * np.array([cos(ccw), sin(ccw)])
                    tri_pts = [sensor_xy, sensor_xy+vec_cw, sensor_xy+vec_ccw]
                    cones.append(plt.Polygon(tri_pts, color=c[i], alpha=0.15))

            # Add this step's artists to the list of artists.
            robots = PatchCollection(robots, match_original=True, zorder=3)
            lines = LineCollection(lines, linewidths=0.5, colors=c, alpha=0.75,\
                                   zorder=2)
            cones = PatchCollection(cones, match_original=True, zorder=1)
            ims.append([title, ax.add_collection(robots), \
                        ax.add_collection(lines), ax.add_collection(cones)])

        # Animate.
        ani = ArtistAnimation(fig, ims, interval=interval, blit=True)
        ani.save('anis/' + self.fname + '_ani' + anno + '.mp4')
        plt.close()


def load_exp(fname):
    """
    Load an experiment from the specified file.
    """
    with open(fname, 'rb') as f:
        exp = pickle.load(f)

    return exp


### DATA EXPERIMENTS ###

def exp_base(seed=None):
    """
    With default parameters, investigate aggregation over time.
    """
    params = {}  # This uses all default values.
    exp = Experiment('base', params, seed=seed)
    exp.run()
    exp.save()
    exp.plot_evo(runs=[0], iters=[0])
    exp.animate(run=0, iter=0)


def exp_symm(seed=None):
    """
    With default parameters and symmetric initialization, investigate
    aggregation over time for a few system sizes.
    """
    N = [3, 5, 10]
    params = {'N' : N, 'init' : ['symm']}
    exp = Experiment('symm', params, seed=seed)
    exp.run()
    exp.save()
    exp.plot_evo(runs=np.arange(len(exp.params)), iters=[0], metrics=['disp'], \
                 labels=['{} robots'.format(i) for i in N], \
                 title='Symmetric Initial Configuration')


def exp_errprob(seed=None):
    """
    With default parameters and a range of error probabilities, investigate
    average time to aggregation with a 15% stopping condition.
    """
    N = [10, 25, 50, 100]
    errprob = np.arange(0, 0.501, 0.0125)
    params = {'N' : N, 'noise' : [('err', p) for p in errprob], 'stop' : [0.15]}
    exp = Experiment('errprob', params, iters=25, savehist=False, seed=seed)
    exp.run()
    exp.save()
    exp.plot_aggtime(N, errprob, 'Error Probability')


def exp_motion(seed=None):
    """
    With default parameters and a range of motion noise strengths, investigate
    average time to aggregation with a 15% stopping condition.
    """
    N = [10, 25, 50, 100]
    fmax = np.arange(0, 40.1, 1.25)
    params = {'N' : N, 'noise' : [('mot', f) for f in fmax], 'stop' : [0.15]}
    exp = Experiment('motion', params, iters=25, savehist=False, seed=seed)
    exp.run()
    exp.save()
    exp.plot_aggtime(N, fmax, 'Max. Noise Force (N)')


def exp_cone(seed=None):
    """
    With default parameters and a range of sight sensor sizes, investigate
    average time to aggregation with a 15% stopping condition.
    """
    N = [10, 25, 50, 100]
    sensor = np.arange(0, np.pi, 0.1)
    params = {'N' : N, 'sensor' : sensor, 'stop' : [0.15]}
    exp = Experiment('cone', params, iters=25, savehist=False, seed=seed)
    exp.run()
    exp.save()
    exp.plot_aggtime(N, sensor, 'Sight Sensor Size (rad)')


### CALIBRATION EXPERIMENTS ###

def exp_step(seed=None):
    """
    With default parameters and a range of time step durations, investigate
    aggregation over time.
    """
    step = [0.0005, 0.001, 0.005, 0.01, 0.025]
    params = {'N' : [50], 'time' : [120], 'step' : step}
    exp = Experiment('step', params, seed=seed)
    exp.run()
    exp.save()
    exp.plot_evo(runs=np.arange(len(exp.params)), iters=[0], metrics=['disp'], \
                 labels=['{}s'.format(i) for i in step])


if __name__ == '__main__':
    # Parse command line arguments.
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-E', '--exps', type=str, nargs='+', required=True, \
                        help='IDs of experiments to run')
    parser.add_argument('-R', '--rand_seed', type=int, default=None, \
                        help='Seed for random number generation')
    args = parser.parse_args()

    # Run selected experiments.
    exps = {'base' : exp_base, 'symm' : exp_symm, 'errprob' : exp_errprob, \
            'motion' : exp_motion, 'cone' : exp_cone, 'step' : exp_step}
    for id in args.exps:
        exps[id](args.rand_seed)
