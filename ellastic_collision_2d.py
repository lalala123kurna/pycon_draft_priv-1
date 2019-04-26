import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anm                                       
import math as mt
import pdb


def advection_1d(x_i, v, dt=1):
    return x_i + v * dt


def collision_1d(v1_i, v2_i, m1=1, m2=1):
    v1_f = (m1 - m2) / (m1 + m2) * v1_i + 2 * m2 / (m1 + m2) * v2_i
    v2_f = (m2 - m1) / (m1 + m2) * v2_i + 2 * m1 / (m1 + m2) * v1_i
    return v1_f, v2_f


def collision_2d(v1, v2, dist, m1, m2):
    n_vect = dist / np.linalg.norm(dist)

    print(dist)
    print(n_vect)

    # scalars
    v1n = np.dot(v1, n_vect)
    v2n = np.dot(v2, n_vect)
    v1s = v1 - v1n * n_vect
    v2s = v2 - v2n * n_vect

    v1n_f, v2n_f = collision_1d(v1n, v2n, m1, m2)
    v1n_f *= n_vect
    v2n_f *= n_vect

    v1_f = v1s + v1n_f
    v2_f = v2s + v2n_f
    # pdb.set_trace()
    return v1_f, v2_f


def simulation_step(dt, m1, m2, radius, loc1_0, loc2_0, v1_0, v2_0, domain):
    """ returns positions and velocity  of the particles after one step"""

    # radius
    r1 = r2 = radius

    # positions  = (x, y)
    loc1 = loc1_0
    loc2 = loc2_0
    dist = loc2 - loc1

    # velocities
    v1 = v1_0
    v2 = v2_0

    # collisions
    if np.linalg.norm(dist) < (r1 + r2):
        v1, v2 = collision_2d(v1, v2, dist, m1, m2)

    # advection
    for i in [0, 1]:
        loc1[i] = advection_1d(x_i=loc1[i], v=v1[i], dt=dt)
        loc2[i] = advection_1d(x_i=loc2[i], v=v2[i], dt=dt)

    # boundary condition
    for i in [0, 1]:
        if loc1[i] - r1 < domain[i][0]:
            v1[i] *= -1
        if loc1[i] + r1 > domain[i][1]:
            v1[i] *= -1
        if loc2[i] - r2 < domain[i][0]:
            v2[i] *= -1
        if loc2[i] + r2 > domain[i][1]:
            v2[i] *= -1

    return [loc1, loc2], [v1, v2]


class Movie_2d:
    """ TODO """

    def __init__(self, dt, t_max, loc1, loc2, vel1, vel2, domain, m1=1, m2=1, radius=10):
        self.dt = dt
        self.fps = int(mt.floor(1/dt))
        self.n_max = int(mt.floor(t_max / dt))
        self.loc1 = loc1                                                         
        self.vel1 = vel1
        self.loc2 = loc2                                                         
        self.vel2 = vel2
        self.domain = domain
        self.m1 = m1
        self.m2 = m2
        self.radius = radius

        self.fig = plt.figure()                                                
        self.fig.subplots_adjust(left=0, bottom=0.1, right=1, top=.9)             

        self.ax = self.fig.add_subplot(111, aspect='equal', autoscale_on=False, xlim=(domain[0]), ylim=domain[1])
        self.c1 = plt.Circle(self.loc1, self.radius, color='r', fill=False, linewidth=2)
        self.c2 = plt.Circle(self.loc2, self.radius, color='r', fill=False, linewidth=2)
        self.ax.add_patch(self.c1)
        self.ax.add_patch(self.c2)

        self.frame = plt.Rectangle([0,0], domain[0][1] - domain[0][0], domain[1][1] - domain[1][0], fc='none')
        self.ax.add_patch(self.frame)        

    def generate_frame(self, i):                                               
        loc_arr, vel_arr = simulation_step(self.dt, self.m1, self.m2, self.radius,
                                           self.loc1, self.loc2, self.vel1, self.vel2, self.domain)
        self.loc1 = loc_arr[0]
        self.loc2 = loc_arr[1]
        self.vel1 = vel_arr[0]
        self.vel2 = vel_arr[1]
        self.c1.set_center(self.loc1)
        self.c2.set_center(self.loc2)
        # TODO: not sure what should be instead of self.c1 and why
        return self.c1, self.c2, self.frame
                                                                               
    def animate(self, name):                                                   
        animation = anm.FuncAnimation(self.fig, self.generate_frame, frames=self.n_max, interval=10, blit=True)
                                                                               
        fps = self.fps                                                         
        animation.save(name+'.mp4', fps=fps, extra_args=['-vcodec', 'libx264'])


# initial values and constants
loc1_0 = np.array([30., 50.])
loc2_0 = np.array([30 + (1200) ** 0.5, 50. + 20])
# loc2_0 = np.array([30 + (1200) ** 0.5, 50.])
v1 = [12., 12.]
v2 = [-8, -15.]
# v1 = [12., 0]
# v2 = [-8, 0]
domain = [[0, 100], [0, 100]]
m1 = 1
m2 = 2
dt = 1/30.
t_max = 1#50
radius = 10

movie = Movie_2d(dt=dt, t_max=t_max, loc1=loc1_0, loc2=loc2_0,
                 vel1=v1, vel2=v2, domain=domain, m1=m1, m2=m2, radius=radius)
movie.animate("test_2d")   
