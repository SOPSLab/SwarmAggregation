# SwarmAggregation

SwarmAggregation is a lightweight, flexible Python simulator for the [Gauci et al. (2014)](https://doi.org/10.1177/0278364914525244) aggregation algorithm for e-puck robot swarms.
Each robot has a sight sensor that returns exactly one bit of information: whether or not another robot is seen.
Informally, the algorithm has two rules:

1. If another robot is seen, then rotate clockwise in place.
2. Otherwise, rotate clockwise about a center of rotation that is 90 degrees counter-clockwise from the sight sensor's axis.

SwarmAggregation is used in [this whitepaper](arXiv TODO) to investigate the robustness of the original algorithm to different types of noise and error as well as its performance when using different types of sight sensors.
Its implementation has the following parameters:

| Parameter | Description |
| --- | --- |
| `N` | Number of robots |
| `R` | Radius of rotation about the center of rotation (m) |
| `r` | Radius of a robot (m) |
| `m` | Mass of a robot (kg) |
| `w0` | Rotational speed of a robot about its center of rotation (rad/s) |
| `w1` | Rotational speed of a robot in place (rad/s) |
| `sensor` | Size of the sight sensor (rad) |
| `noise` | `('err', p)` for sensor error probability with probability `p` |
|  | `('mot', f)` for motion noise with maximum force `f` (N) |
| `time` | Wall-clock duration of a simulation (s) |
| `step` | Wall-clock duration of a time step (s) |


## Getting Started

1. You'll need a command line (Unix-based, Windows Command Prompt, or macOS Terminal) and any Python installation version 3.5 or newer. You will also need the [numpy](https://numpy.org/install/), [scipy](https://www.scipy.org/install.html), [matplotlib](https://matplotlib.org/stable/users/installing.html), and [tqdm](https://github.com/tqdm/tqdm#installation) packages.

2. Clone this repository or download the latest [release](https://github.com/SOPSLab/SwarmAggregation/releases).

3. Create `data/`, `figs/`, and `anis/` directories in the code directory.

4. To reproduce the data and figures from our paper, run the experiments script:
```
python exp.py -E <id_of_experiment> -R <random_seed>
```
A list of experiments and their IDs can be found in `exp.py`. All results from the paper were obtained with seed `3121127542`.

5. While less convenient, you can alternatively run a single simulation of the ARM with your own parameters of interest:
```
python
>>> from aggregation import aggregation
>>> # set up parameters as listed above
>>> # returns the history of configurations and the final step index
>>> history, final = aggregation(N, R, r, m, w0, w1, sensor, noise, time, step)
```

6. Note that the plotting and analysis functions are members of the `Experiment` class in `exp.py`. If you want to do anything heavier than simply getting the data for a single run, you should add your own experiments to `exp.py` and run them as in Step 4.


## Contributing

If you'd like to leave feedback, feel free to open a [new issue](https://github.com/SOPSLab/SwarmAggregation/issues/new/).
If you'd like to contribute, please submit your code via a [pull request](https://github.com/SOPSLab/SwarmAggregation/pulls).
