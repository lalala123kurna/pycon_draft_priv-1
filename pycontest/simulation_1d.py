import numpy
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import matplotlib.animation as anm                                       
import math as mt
import pdb

def advection_1d(x_i, v, dt):
    return x_i + v * dt


def collision_1d(v1_i, v2_i, m1=1, m2=1):
    v1_f = (m1 - m2) / (m1 + m2) * v1_i + 2 * m2 / (m1 + m2) * v2_i
    v2_f = (m2 - m1) / (m1 + m2) * v2_i + 2 * m1 / (m1 + m2) * v1_i
    return v1_f, v2_f


def simulation_step(dt, x1_0, x2_0, v1_0, v2_0, domain_x, radius):
    """ returns positions and velocity  of the particles after one step"""

    # radius
    r1 = r2 = radius

    # positions
    x1 = x1_0
    x2 = x2_0

    # velocities
    v1 = v1_0
    v2 = v2_0

    # do collision between 2 particles
    x1 = advection_1d(x1, dt, v1)
    x2 = advection_1d(x2, dt, v2)
    if abs(x1 - x2) < (r1 + r2):
        v1, v2 = collision_1d(v1, v2)

    #boundary condition        
    if x1 - r1 < domain_x[0]:
        v1 *= -1
    if x1 + r1 > domain_x[1]:
        v1 *= -1
    if x2 - r2 < domain_x[0]:
        v2 *= -1
    if x2 + r2 > domain_x[1]:
        v2 *= -1

    return [x1, x2], [v1, v2]


class Movie:
    """ Movie class

        __init__ input:
            - dt       : simulation timestep
            - t_max    : maximum simulation time
            - loc      : array with balls positions
            - vel      : array with balls velocities
            - domain_x : domain extent

        methods:
            - generate frame() : create a single animation frame
            - animate(name)    : create the animation; name is the string with file name for the animation
    """

    def __init__(self, dt, t_max, loc, vel, domain_x, radius=1):

        self.dt = dt
        self.fps = int(mt.floor(1/dt))
        self.n_max = int(mt.floor(t_max / dt))
        self.loc = loc
        self.vel = vel
        self.domain_x = domain_x
        self.radius = radius

        self.fig = plt.figure()
        self.fig.subplots_adjust(left=0.1, bottom=0, right=.9, top=1)

        self.ax = self.fig.add_subplot(111, aspect='equal', autoscale_on=False,xlim=(domain_x), ylim=(0,3))
        self.c1 = plt.Circle((self.loc[0], 1.5), self.radius, color='r', fill=False, linewidth=2)
        self.c2 = plt.Circle((self.loc[1], 1.5), self.radius, color='r', fill=False, linewidth=2)
        self.ax.add_patch(self.c1)
        self.ax.add_patch(self.c2)

        self.frame = plt.Rectangle([0,0], domain_x[1] - domain_x[0], 3., fc='none')
        self.ax.add_patch(self.frame)

    def generate_frame(self, i):
        self.loc, self.vel = simulation_step(self.dt, self.loc[0], self.loc[1], self.vel[0], self.vel[1],
                                             self.domain_x, self.radius)
        self.c1.set_center((self.loc[0], 1.5))
        self.c2.set_center((self.loc[1], 1.5))
        return self.c1, self.c2, self.frame

    def animate(self, name):
        animation = anm.FuncAnimation(self.fig, self.generate_frame, frames=self.n_max, interval=10, blit=True)

        fps = self.fps
        animation.save(name+'.mp4', fps=fps, extra_args=['-vcodec', 'libx264'])

movie = Movie(1/30., 20, [1, 5], [1, -1], [0, 10])
movie.animate("test_1d")
