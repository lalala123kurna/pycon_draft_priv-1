import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anm                                       
import math as mt
import pdb
from scipy.spatial.distance import pdist, squareform

def advection_1d(x_i, v, dt=1):
    return x_i + v * dt


def collision_1d(v1_i, v2_i, m1=1, m2=1):
    v1_f = (m1 - m2) / (m1 + m2) * v1_i + 2 * m2 / (m1 + m2) * v2_i
    v2_f = (m2 - m1) / (m1 + m2) * v2_i + 2 * m1 / (m1 + m2) * v1_i
    return v1_f, v2_f


def collision_2d(v1, v2, dist, m1, m2):
    n_vect = dist / np.linalg.norm(dist)

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


def simulation_step(dt, mass, radius, loc, vel, domain):
    """ returns positions and velocities of particles after one step"""

    # TODO - can we do better than this?
    loc_x = np.copy(loc)
    loc_y = np.copy(loc)
    loc_x[:,1] = 1
    loc_y[:,0] = 1

    # find pairs of particles undergoing a collision
    dist   = squareform(pdist(loc))
    dist_x = squareform(pdist(loc_x))
    dist_y = squareform(pdist(loc_y))

    ind1, ind2 = np.where(dist < 2 * radius)
    unique = (ind1 < ind2)
    ind1 = ind1[unique]
    ind2 = ind2[unique]

    # collisions
    for id1, id2 in zip(ind1, ind2):
        vel[id1], vel[id2] = collision_2d(vel[id1], vel[id2], [dist_x[id1,id2], dist_y[id1,id2]], mass[id1], mass[id2])

    # advection
    loc[:,0] = advection_1d(x_i=loc[:,0], v=vel[:,0], dt=dt)
    loc[:,1] = advection_1d(x_i=loc[:,1], v=vel[:,1], dt=dt)

    # find outside domain points ...
    out_left   = (loc[:, 0] < domain[0][0] + radius)
    out_right  = (loc[:, 0] > domain[0][1] - radius)
    out_top    = (loc[:, 1] < domain[1][0] + radius)
    out_bottom = (loc[:, 1] > domain[1][1] - radius)
    
    # TODO ... fix their positions ...?
    #loc[out_left, 0]   = domain[0][0] + radius
    #loc[out_right, 0]  = domain[0][1] - radius
    #loc[out_top, 1]    = domain[1][0] + radius
    #loc[out_bottom, 1] = domain[1][1] - radius

    # ... change their velocities
    vel[out_left | out_right,  0] *= -1
    vel[out_top  | out_bottom, 1] *= -1

    return loc, vel

class Movie_2d:
    """ TODO """

    def __init__(self, dt, t_max, loc, vel, domain, mass, radius=10):
        self.dt = dt
        self.fps = int(mt.floor(1/dt))
        self.n_max = int(mt.floor(t_max / dt))
        self.loc = loc                                                         
        self.vel = vel
        self.domain = domain
        self.mass = mass
        self.radius = radius

        self.fig = plt.figure()                                                
        self.fig.subplots_adjust(left=0, bottom=0.1, right=1, top=.9)             

        self.ax = self.fig.add_subplot(111, aspect='equal', autoscale_on=False, xlim=(domain[0]), ylim=domain[1])
        for el in range(loc.shape[0]):
            self.circ = plt.Circle(self.loc[el,:], self.radius, color='r', fill=False, linewidth=2)
            self.ax.add_patch(self.circ)

        self.frame = plt.Rectangle([0,0], domain[0][1] - domain[0][0], domain[1][1] - domain[1][0], fc='none')
        self.ax.add_patch(self.frame)        

    def generate_frame(self, i):                                               
 
        loc, vel = simulation_step(self.dt, self.mass, self.radius,
                                   self.loc, self.vel, self.domain)
        for el in range(loc.shape[0]):
            self.circ.set_center(self.loc[el,:])
        return self.circ, self.frame

    def animate(self, name):                                                   
        animation = anm.FuncAnimation(self.fig, self.generate_frame, frames=self.n_max, interval=10, blit=True)

        fps = self.fps                                                         
        animation.save(name+'.mp4', fps=fps, extra_args=['-vcodec', 'libx264'])


# initial values and constants
vel = np.array([[12., 12.], [-8, -15.]])
loc = np.array([[30., 50.], [30 + (1200) ** 0.5, 50. + 20]])
domain = np.array([[0., 100.], [0., 100.]])
mass = np.array([1., 2.])
dt = 1/30.
t_max = 1#50
radius = 10
 
movie = Movie_2d(dt=dt, t_max=t_max, loc=loc,
                 vel=vel, domain=domain, mass=mass, radius=radius)
movie.animate("test_2d_2_balls")   


# initial values and constants
vel = np.array([[12., 12.], [-8, -15.], [-3, 4]])
loc = np.array([[30., 50.], [30 + (1200) ** 0.5, 50. + 20], [75, 75]])
domain = np.array([[0., 100.], [0., 100.]])
mass = np.array([1., 2., 1.])
dt = 1/30.
t_max = 50
radius = 10
 
movie = Movie_2d(dt=dt, t_max=t_max, loc=loc,
                 vel=vel, domain=domain, mass=mass, radius=radius)
movie.animate("test_2d_3_balls")   

