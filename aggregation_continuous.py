import datetime
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





def truncate(number, decimals=0):
    if decimals == 0:
        return math.trunc(number)
    factor = 10.0 ** decimals
    return math.trunc(number * factor) / factor




def dist(a, b) :
    dist = math.sqrt( (b.X - a.X)**2 + (b.Y - a.Y)**2 )
    return dist




def init_sys_symmetric(n) :
    robot_array = []
    circ_traj_rot_radius = 14.45
    puck_radius = 3.7
    # initial_config_circ_radius = (n+1) * circ_traj_rot_radius * 3
    initial_config_circ_radius = (n+1) * circ_traj_rot_radius
    if (n==10):
        initial_config_circ_radius = 125
    # initial_config_circ_radius = 45
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




def init_sys_custom():
    # Deadlock - rigid (most basic): one pair of paritcles
    # See data/expdeadlock_seed_datetime_run0_iter0animation.pkl for direct visualization
    coords_thetas = [ [0,0, math.pi], [0,7.4, 0] ]

    # # Deadlock - rigid: one pair of paritcles and one triplet of particles (mirrors figure in paper)
    # # See data/expdeadlock_45907581_190324472155_run0_iter0animation.pkl for direct visualization
    # coords_thetas = [ [0,0, math.pi], [0,7.4, 0], [40,0, 5*math.pi/6], [40,7.4, math.pi/6], [40+(math.sqrt(7.4**2-3.7**2)),3.7, 3*math.pi/2] ]

    # # Deadlock - "natural": two individually aggregated yet separated clusters
    # # See data/expdeadlock_987342234_203227346276_run0_iter0animation.pkl for direct visualization
    # coords_thetas = [ [0,0, math.pi/2], [7.4,0, math.pi], [14.8,0, 0], [22.2,0, 3*math.pi/2],   [3.7,math.sqrt(7.4**2-3.7**2), math.pi/4], [11.1,math.sqrt(7.4**2-3.7**2), 0], [18.5,math.sqrt(7.4**2-3.7**2), 7*math.pi/4],   [3.7,-1*math.sqrt(7.4**2-3.7**2), 3*math.pi/4], [11.1,-1*math.sqrt(7.4**2-3.7**2), math.pi], [18.5,-1*math.sqrt(7.4**2-3.7**2), 5*math.pi/4],
    #                   [100,-35, math.pi/2], [100+7.4,-35, math.pi], [100+14.8,-35, 3*math.pi/2],   [100+3.7,-35+math.sqrt(7.4**2-3.7**2), math.pi/6], [100+11.1,-35+math.sqrt(7.4**2-3.7**2), 11*math.pi/6],   [100+3.7,-35-math.sqrt(7.4**2-3.7**2), 5*math.pi/6], [100+11.1,-35-math.sqrt(7.4**2-3.7**2), 7*math.pi/6] ]

    robot_array = []
    for i in range(len(coords_thetas)):
        robot_array.append( Robot(coords_thetas[i][0], coords_thetas[i][1], coords_thetas[i][2], 3.7, .152, 0, 0) )

    return robot_array










def collision_control(arr, time_step) :
    # spring_constant = 3090
    # spring_constant = 100
    # spring_constant = 16.5

    spring_constant = 21964 / (28.9 * (time_step**2))

    # overlap_tmp = 7.4 - math.sqrt( ( 7.4-28.9*math.sin(.00075*time_step) )**2 + ( 28.9*(math.cos(-.00075*time_step)-1) )**2 )
    # theta_tmp = math.atan( (-7.4-28.9*math.sin(-.00075*time_step)) / (28.9-28.9*math.cos(-.00075*time_step)) )
    # spring_constant = abs( (21964*math.sin(.00075*time_step)) / (time_step**2 * math.sin(theta_tmp) * overlap_tmp) )

    # spring_constant_2 = (21964*(math.cos(math.pi-.00075*time_step)+1)) / (time_step**2 * math.cos(theta_tmp) * overlap_tmp)
    # print(spring_constant)
    # print(spring_constant_2)

    # spring_constant = (21964*math.sin(math.pi/4-.00075*time_step)-21964) / (time_step**2 * math.sin(theta_tmp) * overlap_tmp)
    # print(spring_constant)



    for a in range(len(arr)):
        for b in range(a+1, len(arr)):
            if(dist(arr[a], arr[b]) > 0 and dist(arr[a], arr[b]) < 7.4):
                spring_force = spring_constant * ( (2*arr[a].robot_radius) - (dist(arr[a], arr[b])) )

                if (arr[b].X-arr[a].X == 0):
                    if (arr[b].Y > arr[a].Y):
                        theta = math.pi/2
                    else:
                        theta = 3*math.pi/2
                else:
                    theta = math.atan( (arr[b].Y-arr[a].Y) / (arr[b].X-arr[a].X) )

                x_compon_force = spring_force * math.cos(theta)
                y_compon_force = spring_force * math.sin(theta)

                if (arr[b].X > arr[a].X and arr[b].Y < arr[a].Y):  # -x,+y
                    if (x_compon_force > 0):
                        x_compon_force = x_compon_force * (-1)
                    if (y_compon_force < 0):
                        y_compon_force = y_compon_force * (-1)
                elif (arr[b].X < arr[a].X and arr[b].Y < arr[a].Y):   # +x,+y
                    if (x_compon_force < 0):
                        x_compon_force = x_compon_force * (-1)
                    if (y_compon_force < 0):
                        y_compon_force = y_compon_force * (-1)
                elif (arr[b].X > arr[a].X and arr[b].Y > arr[a].Y):   # -x,-y
                    if (x_compon_force > 0):
                        x_compon_force = x_compon_force * (-1)
                    if (y_compon_force > 0):
                        y_compon_force = y_compon_force * (-1)
                elif (arr[b].X < arr[a].X and arr[b].Y > arr[a].Y):   # +x,-y
                    if (x_compon_force < 0):
                        x_compon_force = x_compon_force * (-1)
                    if (y_compon_force > 0):
                        y_compon_force = y_compon_force * (-1)

                arr[a].forces_x += x_compon_force
                arr[a].forces_y += y_compon_force
                arr[b].forces_x += (-1)*x_compon_force
                arr[b].forces_y += (-1)*y_compon_force










def motion_noise(robot_array, max_force, seed):
    rng = np.random.default_rng(seed)

    for i in range(len(robot_array)):
        motion_force = max_force * rng.random()
        direction = 2*math.pi * rng.random()
        robot_array[i].forces_x += motion_force * math.cos(direction)
        robot_array[i].forces_y += motion_force * math.sin(direction)





def elastic_poten_and_motion_to_kinetic(robot_array, _time_step):
    for i in range(len(robot_array)):
        x_impulse = robot_array[i].forces_x * (_time_step/1000)      # N*s or kg*m/s
        x_velo = x_impulse / robot_array[i].mass        # m/s
        converted_cm_ms_x_velo = x_velo / 10
        robot_array[i].X += converted_cm_ms_x_velo * _time_step

        y_impulse = robot_array[i].forces_y * (_time_step/1000)     # N*s or kg*m/s
        y_velo = y_impulse / robot_array[i].mass        # m/s
        converted_cm_ms_y_velo = y_velo / 10
        robot_array[i].Y += converted_cm_ms_y_velo * _time_step




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








def update_system(prev_step_arr, errorprob_bool, noise_val, sensor_alpha, _time_step, seed) :
    rng = np.random.default_rng(seed)

    prev_system_arr = prev_step_arr
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
            new_theta = (prev_system_arr[i].theta - (.00502*_time_step)) % (2*math.pi)
            if (new_theta < 0) :
                new_theta = (2*math.pi) + new_theta
            new_pos_x = prev_system_arr[i].X
            new_pos_y = prev_system_arr[i].Y
        else :      # no robot in sight -> move clockwise on circ. trajectory with R = 14.45 cm @ -.75 rad/s
            new_pos_x = prev_system_arr[i].X + ( circ_traj_radius*math.cos(prev_system_arr[i].theta-(.00075*_time_step)) - circ_traj_radius*math.cos(prev_system_arr[i].theta) )
            new_pos_y = prev_system_arr[i].Y + ( circ_traj_radius*math.sin(prev_system_arr[i].theta-(.00075*_time_step)) - circ_traj_radius*math.sin(prev_system_arr[i].theta) )

            new_theta = (prev_system_arr[i].theta - (.00075*_time_step)) % (2*math.pi)
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

    # print(dispersion_val)

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
        return (True, dispersion_val / ideal_cluster)     # TERMINATE RUN; STOPPING CONDITION MET
    else:
        return (False, dispersion_val / ideal_cluster)









def aggregation(_N, _T, _init, _noise, _time_step, _stopping, _savehistory, _seed):
    # start_time = datetime.datetime.now()

    if _init is 'random':
        init_config = init_sys_random(_N, _seed)
    elif _init is 'symmetric':
        init_config = init_sys_symmetric(_N)
    elif _init is 'custom':
        init_config = init_sys_custom()
    else:
        assert False, 'ERROR: Unrecognized initialization method: ' + _init

    if (_savehistory == False):
        step = 1
        curr_step = init_config
        step += 1
        while (step < 300 * (1/_time_step)*1000):
            if (_noise[0] == 'errorprob'):
                curr_step = update_system(curr_step, True, _noise[1], _T, _time_step, _seed)
            elif (_noise[0] == 'motion'):
                curr_step = update_system(curr_step, False, 0, _T, _time_step, _seed)
                motion_noise(curr_step, _noise[1], _seed)
            else:
                assert False, 'ERROR: Unrecognized interaction rule ' + _noise[0]
            collision_control(curr_step, _time_step)
            elastic_poten_and_motion_to_kinetic(curr_step, _time_step)
            clear_forces(curr_step)
            if (step % (.5 * 1000*(1/_time_step)) == 0):
                stopping = dispersion_stopping_condition(curr_step)
                if (stopping[0] == True):
                    break
            step += 1
        # print(step)
        # print(stopping[1])
        # print(datetime.datetime.now() - start_time)
        return step
    else :   # _savehistory == True
        if (_stopping[0] == True):
            history = []
            history.append(init_config)
            step = 1
            while (step < 300 * (1/_time_step)*1000):
                if (_noise[0] == 'errorprob'):
                    history.append( update_system(history[step-1], True, _noise[1], _T, _time_step, _seed) )
                elif (_noise[0] == 'motion'):
                    history.append( update_system(history[step-1], False, 0, _T, _time_step, _seed) )
                    motion_noise(history[step], _noise[1], _seed)
                else:
                    assert False, 'ERROR: Unrecognized interaction rule ' + _noise[0]
                collision_control(history[step], _time_step)
                elastic_poten_and_motion_to_kinetic(history[step], _time_step)
                clear_forces(history[step])
                if (step % (1 * 1000*(1/_time_step)) == 0):
                    stopping = dispersion_stopping_condition(history[step])
                    if (stopping[0] == True):
                        break
                step += 1
            allRawData_arr = robotObjArr_to_RawDataArr(history)
            return np.array(allRawData_arr)
        else:   # _stopping[0] == False
            history = []
            history.append(init_config)
            for step in range(1, _stopping[1]):
                if (_noise[0] == 'errorprob'):
                    history.append( update_system(history[step-1], True, _noise[1], _T, _time_step, _seed) )
                elif (_noise[0] == 'motion'):
                    history.append( update_system(history[step-1], False, 0, _T, _time_step, _seed) )
                    motion_noise(history[step], _noise[1], _seed)
                else:
                    assert False, 'ERROR: Unrecognized interaction rule ' + _noise[0]
                collision_control(history[step], _time_step)
                elastic_poten_and_motion_to_kinetic(history[step], _time_step)
                clear_forces(history[step])
            allRawData_arr = robotObjArr_to_RawDataArr(history)
            # print(datetime.datetime.now() - start_time)
            return np.array(allRawData_arr)
