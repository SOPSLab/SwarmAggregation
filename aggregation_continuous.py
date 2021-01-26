# Project:     Swarm Aggregation

import math
import numpy as np
import random




class Robot:
    def __init__(self, X, Y, theta, robot_radius, mass, forces_x, forces_y) :
        # x position, cm
        self.X = X

        # y position, cm
        self.Y = Y

        # angle in radians of the line segment starting @ center of rotation pointing towards the robot; 0, 2pi rad is due east; [0, 2pi)
        # center of rotation is the "origin" of the "unit circle"
        self.theta = theta

        # radius of the puck robot itself, cm
        self.robot_radius = robot_radius

        # mass in kg or robot puck
        self.mass = mass

        # sum of all forces acting on robot resulting from "elastic collisions" and "motion noise"; x and y components separately
        self.forces_x = forces_x
        self.forces_y = forces_y




def dist(a, b) :
    dist = math.sqrt( (b.X - a.X)**2 + (b.Y - a.Y)**2 )
    return dist




def init_sys_symmetric(n) :
    robot_array = []
    circ_traj_rot_radius = 14.45
    puck_radius = 3.7
    # initial_config_circ_radius = (n+1) * circ_traj_rot_radius * 3
    # initial_config_circ_radius = (n+1) * circ_traj_rot_radius
    initial_config_circ_radius = 45
    # initial_config_circ_radius = 25
    theta_step = (2*math.pi) / n
    for i in range(0, n) :
        theta_around_init_config_circ = theta_step * i
        x_coor = initial_config_circ_radius * math.cos(theta_around_init_config_circ)
        y_coor = initial_config_circ_radius * math.sin(theta_around_init_config_circ)
        line_sensor_theta = (  ( (i+1)*theta_step ) + (.5 * ( ((n-2)*math.pi) / n ) ) ) % (2*math.pi)
        theta = (line_sensor_theta - .5*math.pi) % (2*math.pi)
        robot_array.append(Robot(x_coor, y_coor, theta, puck_radius, .152, 0, 0))
    return robot_array


def init_sys_random(n, seed) :
    robot_array = []
    robot_count = 0
    puck_radius = 3.7
    box_radius = ( math.sqrt( (puck_radius**2 * math.pi * n) / .019 ) ) / 2

    rng = np.random.default_rng(seed)

    while (robot_count < n):
        x_coor = ((rng.random()-.5)*2) * box_radius
        y_coor = ((rng.random()-.5)*2) * box_radius
        theta = rng.random() * 2*math.pi
        robot_tmp = Robot(x_coor, y_coor, theta, puck_radius, .152, 0, 0)

        init_config_valid_bool = True
        for x in range(0, len(robot_array)):
            if (dist(robot_tmp, robot_array[x]) <= 2*puck_radius):
                init_config_valid_bool = False
                break

        if (init_config_valid_bool == True):
            robot_array.append(robot_tmp)
            robot_count += 1

    return robot_array




def collision_control(arr) :
    spring_constant = 6080
    for a in range(len(arr)):
        for b in range(len(arr)):
            robot = arr[a]
            other = arr[b]
            if(dist(robot, other) > 0 and dist(robot, other) < (2*robot.robot_radius)):
                # print("collision.")
                spring_force = spring_constant * ( (2*robot.robot_radius) - (dist(robot, other)) )
                theta = math.atan( (other.Y-robot.Y) / (other.X-robot.X) )
                x_compon_force = spring_force * math.cos(theta)
                y_compon_force = spring_force * math.sin(theta)
                if (other.X > robot.X and other.Y > robot.Y):
                    x_compon_force = x_compon_force * (-1)
                    y_compon_force = y_compon_force * (-1)
                elif (other.X < robot.X and other.Y < robot.Y):
                    x_compon_force = x_compon_force
                    y_compon_force = y_compon_force
                elif (other.X > robot.X and other.Y < robot.Y):
                    x_compon_force = x_compon_force * (-1)
                    y_compon_force = y_compon_force * (-1)
                elif (other.X < robot.X and other.Y > robot.Y):
                    x_compon_force = x_compon_force
                    y_compon_force = y_compon_force
                robot.forces_x += x_compon_force
                robot.forces_y += y_compon_force
        arr[a] = robot




def motion_noise(robot_array, max_force, seed):
    rng = np.random.default_rng(seed)

    for i in range(len(robot_array)):
        motion_force = max_force * rng.random()
        direction = 2*math.pi * rng.random()
        robot_array[i].forces_x += motion_force * math.cos(direction)
        robot_array[i].forces_y += motion_force * math.sin(direction)




def elastic_poten_and_motion_to_kinetic(robot_array):
    for i in range(len(robot_array)):
        x_impulse = robot_array[i].forces_x * .0005      # N*s or kg*m/s
        x_velo = x_impulse / robot_array[i].mass        # m/s
        converted_cm_ms_x_velo = x_velo / 10
        robot_array[i].X += (converted_cm_ms_x_velo / 2)

        y_impulse = robot_array[i].forces_y * .0005     # N*s or kg*m/s
        y_velo = y_impulse / robot_array[i].mass        # m/s
        converted_cm_ms_y_velo = y_velo / 10
        robot_array[i].Y += (converted_cm_ms_y_velo / 2)




def clear_forces(robot_array):
    for i in range(len(robot_array)):
        robot_array[i].forces_x = 0
        robot_array[i].forces_y = 0




def robot_in_sight(robot, robot_array, alpha):
    line_sensor_theta = ( robot.theta + (math.pi/2) ) % (2*math.pi)      # angle in radians of the direction the line sensor in pointed towards; 0, 2pi is due east; [0, 2pi)

    xr = robot.X
    yr = robot.Y

    rr = robot.robot_radius

    for i in range(0, len(robot_array)):
        x2 = robot_array[i].X
        y2 = robot_array[i].Y

        curr_robot_in_sight_bool = True

        if ( line_sensor_theta+(alpha/2) > 3*math.pi/2 or line_sensor_theta+(alpha/2) < math.pi/2 ):
            in_sight_bool_tmp = y2 <= math.tan(line_sensor_theta+(alpha/2))*(x2-xr) + yr + rr/math.cos(line_sensor_theta+(alpha/2))
            if (in_sight_bool_tmp == False):
                curr_robot_in_sight_bool = False
        elif ( line_sensor_theta+(alpha/2) > math.pi/2 and line_sensor_theta+(alpha/2) < 3*math.pi/2 ):
            in_sight_bool_tmp = y2 >= math.tan(line_sensor_theta+(alpha/2))*(x2-xr) + yr + rr/math.cos(line_sensor_theta+(alpha/2))
            if (in_sight_bool_tmp == False):
                curr_robot_in_sight_bool = False
        elif ( line_sensor_theta+(alpha/2) == math.pi/2 ):
            in_sight_bool_tmp = x2 >= xr-rr
            if (in_sight_bool_tmp == False):
                curr_robot_in_sight_bool = False
        else:   # line_sensor_theta+(alpha/2) == 3*math.pi/2
            in_sight_bool_tmp = x2 <= xr+rr
            if (in_sight_bool_tmp == False):
                curr_robot_in_sight_bool = False

        if ( line_sensor_theta-(alpha/2) > 3*math.pi/2 or line_sensor_theta-(alpha/2) < math.pi/2 ):
            in_sight_bool_tmp = y2 >= math.tan(line_sensor_theta-(alpha/2))*(x2-xr) + yr - rr/math.cos(line_sensor_theta-(alpha/2))
            if (in_sight_bool_tmp == False):
                curr_robot_in_sight_bool = False
        elif ( line_sensor_theta-(alpha/2) < 3*math.pi/2 and line_sensor_theta-(alpha/2) > math.pi/2 ):
            in_sight_bool_tmp = y2 <= math.tan(line_sensor_theta-(alpha/2))*(x2-xr) + yr - rr/math.cos(line_sensor_theta-(alpha/2))
            if (in_sight_bool_tmp == False):
                curr_robot_in_sight_bool = False
        elif ( line_sensor_theta-(alpha/2) == math.pi/2 ):
            in_sight_bool_tmp = x2 <= xr+rr
            if (in_sight_bool_tmp == False):
                curr_robot_in_sight_bool = False
        else:   # line_sensor_theta-(alpha/2) == 3*math.pi/2
            in_sight_bool_tmp = x2 >= xr-rr
            if (in_sight_bool_tmp == False):
                curr_robot_in_sight_bool = False

        if ( line_sensor_theta < math.pi and line_sensor_theta != 0 ):
            in_sight_bool_tmp = y2 > math.tan((math.pi/2)+line_sensor_theta)*(x2-xr) + yr
            if (in_sight_bool_tmp == False):
                curr_robot_in_sight_bool = False
        elif ( line_sensor_theta > math.pi and line_sensor_theta < 2*math.pi ):
            in_sight_bool_tmp = y2 < math.tan((math.pi/2)+line_sensor_theta)*(x2-xr) + yr
            if (in_sight_bool_tmp == False):
                curr_robot_in_sight_bool = False
        elif ( line_sensor_theta == 0 or line_sensor_theta == 2*math.pi ):
            in_sight_bool_tmp = x2 > xr
            if (in_sight_bool_tmp == False):
                curr_robot_in_sight_bool = False
        else:   # line_sensor_theta == math.pi
            in_sight_bool_tmp = x2 < xr
            if (in_sight_bool_tmp == False):
                curr_robot_in_sight_bool = False

        if (curr_robot_in_sight_bool == True):
            return True

    return False








def update_system(prev_system_arr, errorprob_bool, noise_val, sensor_alpha, seed) :
    rng = np.random.default_rng(seed)

    new_system_arr = []

    circ_traj_radius = 14.45        # Gauci paper; in cm

    for i in range(0, len(prev_system_arr)) :
        sight_scope = robot_in_sight(prev_system_arr[i], prev_system_arr, sensor_alpha)

        if (errorprob_bool == True):    # error prob. noise
            prob_num = rng.random()
            if (prob_num < noise_val) :
                if (sight_scope == True):
                    sight_scope = False
                else:
                    sight_scope = True

        if (sight_scope == True) :      # robot in sight -> rotate clockwise in place @ -5.02 rad/s
            # print("robot in sight")
            # new_theta = (prev_system_arr[i].theta - 5*(.00502/2)) % (2*math.pi)
            new_theta = (prev_system_arr[i].theta - (.00502/2)) % (2*math.pi)
            if (new_theta < 0) :
                new_theta = (2*math.pi) + new_theta
            new_pos_x = prev_system_arr[i].X
            new_pos_y = prev_system_arr[i].Y
        else :      # no robot in sight -> move clockwise on circ. trajectory with R = 14.45 cm @ -.75 rad/s
            new_pos_x = prev_system_arr[i].X + ( circ_traj_radius*math.cos(prev_system_arr[i].theta-(.00075/2)) - circ_traj_radius*math.cos(prev_system_arr[i].theta) )
            new_pos_y = prev_system_arr[i].Y + ( circ_traj_radius*math.sin(prev_system_arr[i].theta-(.00075/2)) - circ_traj_radius*math.sin(prev_system_arr[i].theta) )

            new_theta = (prev_system_arr[i].theta - (.00075/2)) % (2*math.pi)
            if (new_theta < 0) :
                new_theta = (2*math.pi) + new_theta

        new_system_arr.append(Robot(new_pos_x, new_pos_y, new_theta, 3.7, .152, 0, 0))

    return new_system_arr




def robotObjArr_to_RawDataArr(robotarr):
    allRawData = []
    for step in range(len(robotarr)):
        currStep = []
        for robot in range(len(robotarr[step])):
            currStep.append( [ robotarr[step][robot].X, robotarr[step][robot].Y, robotarr[step][robot].theta ] )
        allRawData.append(currStep)
    return allRawData





def distArrays(a, b):
    return math.sqrt((b[0] - a[0])**2 + (b[1]-a[1])**2)


def dispersion_stopping_condition(robotarr):
    n = len(robotarr)
    rr = robotarr[0].robot_radius

    x_sum = 0
    y_sum = 0
    for i in range(n):
        x_sum += robotarr[i].X
        y_sum += robotarr[i].Y
    centroid = [x_sum / n, y_sum / n]
    dispersion_val = 0
    for k in range(n):
        dispersion_val += distArrays([robotarr[k].X, robotarr[k].Y], centroid)

    if (n >= 169):
        ideal_cluster = (1*(2*rr) * 6) + (2*(2*rr) * 12) + (3*(2*rr) * 18) + (4*(2*rr) * 24) + (5*(2*rr) * 30) + (6*(2*rr) * 36) + (7*(2*rr) * 42) + (8*(2*rr) * (n-169))
    elif (n >= 127):
        ideal_cluster = (1*(2*rr) * 6) + (2*(2*rr) * 12) + (3*(2*rr) * 18) + (4*(2*rr) * 24) + (5*(2*rr) * 30) + (6*(2*rr) * 36) + (7*(2*rr) * (n-127))
    elif (n >= 91):
        ideal_cluster = (1*(2*rr) * 6) + (2*(2*rr) * 12) + (3*(2*rr) * 18) + (4*(2*rr) * 24) + (5*(2*rr) * 30) + (6*(2*rr) * (n-91))
    elif (n >= 61):
        ideal_cluster = (1*(2*rr) * 6) + (2*(2*rr) * 12) + (3*(2*rr) * 18) + (4*(2*rr) * 24) + (5*(2*rr) * (n-61))
    elif (n >= 37):
        ideal_cluster = (1*(2*rr) * 6) + (2*(2*rr) * 12) + (3*(2*rr) * 18) + (4*(2*rr) * (n-37))
    elif (n >= 19):
        ideal_cluster = (1*(2*rr) * 6) + (2*(2*rr) * 12) + (3*(2*rr) * (n-19))
    elif (n >= 7):
        ideal_cluster = (1*(2*rr) * 6) + (2*(2*rr) * (n-7))
    else:
        ideal_cluster = (1*(2*rr) * (n-1))

    if (dispersion_val <= 1.15 * ideal_cluster):
    # if (dispersion_val <= 1.20 * ideal_cluster):
        return (True, dispersion_val / ideal_cluster)     # TERMINATE RUN; STOPPING CONDITION MET
    else:
        return (False, dispersion_val / ideal_cluster)







def aggregation(_N, _T, _init, _noise, _stopping, _savehistory, _seed):
    if _init is 'random':
        config = np.array(init_sys_random(_N, _seed))
    elif _init is 'symmetric':
        config = np.array(init_sys_symmetric(_N))
    else:
        assert False, 'ERROR: Unrecognized initialization method: ' + _init
    init_config = np.copy(config)

    history = []
    history.append(init_config)

    if (_stopping[0] == True):
        step = 1
        while (True):
            if _noise[0] == 'errorprob':
                history.append( update_system(history[step-1], 1, _noise[1], _T, _seed) )
                collision_control(history[step])
                elastic_poten_and_motion_to_kinetic(history[step])
                clear_forces(history[step])
            elif _noise[0] == 'motion':
                history.append( update_system(history[step-1], 0, 0, _T, _seed) )
                collision_control(history[step])
                motion_noise(history[step], _noise[1], _seed)
                elastic_poten_and_motion_to_kinetic(history[step])
                clear_forces(history[step])
            else:
                assert False, 'ERROR: Unrecognized interaction rule ' + _noise[0]

            if (step % 1000*2 == 0):
                stopping = dispersion_stopping_condition(history[step-1])
                if (stopping[0] == False):
                    break
                elif (step % 60 * 1000*2 == 0):
                    print(stopping[1])

            step += 1

    else:
        for step in range(1, _stopping[1]):
            if _noise[0] == 'errorprob':
                history.append( update_system(history[step-1], 1, _noise[1], _T, _seed) )
                # collision_control(history[step])
                elastic_poten_and_motion_to_kinetic(history[step])
                clear_forces(history[step])
            elif _noise[0] == 'motion':
                history.append( update_system(history[step-1], 0, 0, _T, _seed) )
                # collision_control(history[step])
                motion_noise(history[step], _noise[1], _seed)
                elastic_poten_and_motion_to_kinetic(history[step])
                clear_forces(history[step])
            else:
                assert False, 'ERROR: Unrecognized interaction rule ' + _noise[0]


    allRawData_arr = robotObjArr_to_RawDataArr(history)
    if (_savehistory == True):
        return np.array(allRawData_arr)
    else:
        return len(allRawData_arr)
