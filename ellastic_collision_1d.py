import numpy


def position(x_i, v, dt=1):
    return x_i + v * dt


def collision_1d(v1_i, v2_i, m1=1, m2=1):
    v1_f = (m1 - m2) / (m1 + m2) * v1_i + 2 * m2 / (m1 + m2) * v2_i
    v2_f = (m2 - m1) / (m1 + m2) * v2_i + 2 * m1 / (m1 + m2) * v1_i
    return v1_f, v2_f


def simulation_time(nt, dt=1, x1_0=1, x2_0=5, v1_0=1, v2_0=-1):
    """ returns position of the particles after nt steps"""
    x1 = x1_0
    x2 = x2_0
    v1 = v1_0
    v2 = v2_0
    for it in range(nt):
        x1_i, x2_i = x1, x2
        x1 = position(x1_i, dt, v1)
        x2 = position(x2_i, dt, v2)
        if abs(x1-x2) < abs(v1-v2)*dt:
            v1, v2 = collision_1d(v1, v2)
    return x1, x2


def simulation_collision(x1_0, x2_0, v1, v2, dt=1, m1=1, m2=1):
    """returns position and time when particles collide and the final velocities"""
    x1 = x1_0
    x2 = x2_0
    del_x = 1#abs(v1-v2) * dt
    time = 0
    while abs(x1-x2) > del_x:
        x1_i, x2_i = x1, x2
        x1 = position(x1_i, dt, v1)
        x2 = position(x2_i, dt, v2)
        time += dt

    # final positions and time
    x1_f = x1
    x2_f = x2
    time_f = time

    # final velocities after collision
    v1_f, v2_f = collision_1d(v1, v2, m1, m2)
    return x1_f, x2_f, time_f, v1_f, v2_f

