import numpy as np
#import matplotlib
#matplotlib.use("TkAgg")
#import matplotlib.pyplot as plt
#import matplotlib.animation as anm                                       
import math as mt
#import pdb
from scipy.spatial.distance import pdist, squareform

from pycontest.collisions import collision_2d
import pycontest.movies as mv


def simulation_step(dt, mass, radius, loc, vel, domain):
    """ returns positions and velocities of particles after one step"""
    # TODO: should we have different raidus?

    r12_x = squareform(pdist(loc, lambda r1, r2: r1[0] - r2[0]))
    r12_y = squareform(pdist(loc, lambda r1, r2: r1[1] - r2[1]))
    r12 = np.stack((r12_x, r12_y), axis=-1)
    dist = np.apply_along_axis(np.linalg.norm, 2, r12)

    # finding particles that are colliding
    ind1, ind2 = np.where(dist < 2 * radius)
    unique = (ind1 < ind2)

    ind1 = ind1[unique]
    ind2 = ind2[unique]

    # collisions
    for id1, id2 in zip(ind1, ind2):
        vel[id1], vel[id2] = collision_2d(vel[id1], vel[id2], r12[id1, id2], mass[id1], mass[id2])

    # advection
    loc[:] = loc[:] + vel[:] * dt

    # find outside domain points ...
    out_left   = (loc[:, 0] < domain[0][0] + radius)
    out_right  = (loc[:, 0] > domain[0][1] - radius)
    out_top    = (loc[:, 1] < domain[1][0] + radius)
    out_bottom = (loc[:, 1] > domain[1][1] - radius)
    
    # ... change their velocities
    vel[out_left | out_right,  0] *= -1
    vel[out_top  | out_bottom, 1] *= -1

    return loc, vel


if __name__ == '__main__':

    # initial values and constants
    vel = np.array([[12., 12.], [-8, -15.]])
    loc = np.array([[30., 50.], [30 + (1200) ** 0.5, 50. + 20]])
    domain = np.array([[0., 100.], [0., 100.]])
    mass = np.array([1., 2.])
    dt = 1/30.
    t_max = 30
    radius = 10

    movie = mv.Movie_2d(fun=simulation_step, dt=dt, t_max=t_max, loc=loc,
                        vel=vel, domain=domain, mass=mass, radius=radius)
    movie.animate("test_2d_2_balls_new")


    vel = np.array([[12., 12.], [-8, -15.], [-3, 4]])
    loc = np.array([[30., 50.], [30 + (1200) ** 0.5, 50. + 20], [85, 85]])
    domain = np.array([[0., 100.], [0., 100.]])
    mass = np.array([1., 2., 1.])
    dt = 1/30.
    t_max = 30
    radius = 10

    movie = mv.Movie_2d(fun=simulation_step, dt=dt, t_max=t_max, loc=loc,
                        vel=vel, domain=domain, mass=mass, radius=radius)
    movie.animate("test_2d_3_balls_new")


    vel = np.array([[12., 12.], [-8, -15.], [-3, 4], [10, 0]])
    loc = np.array([[30., 50.], [30 + (1200) ** 0.5, 50. + 20], [85, 85], [10,10]])
    domain = np.array([[0., 100.], [0., 100.]])
    mass = np.array([1., 2., 1., 3])
    dt = 1 / 30.
    t_max = 30
    radius = 10

    movie = mv.Movie_2d(fun=simulation_step, dt=dt, t_max=t_max, loc=loc,
                        vel=vel, domain=domain, mass=mass, radius=radius)
    movie.animate("test_2d_4_balls_new")


    vel_many = np.array([[12., 12.], [-8, -15.], [-3, 4], [10, 0], [2., 12.], [-8, -15.],
                       [3, -4], [-10, 0], [-12., 12.], [8, -15.], [-3, -4], [-10, 0]])
    loc_many = np.array([[20., 45.], [65, 70.], [85., 90.], [10., 10.], [50., 45.], [15, 70.],
                    [15., 90.], [45., 10.], [75., 45.], [40, 70.], [45., 90.], [90., 10.]])
    domain = np.array([[0., 100.], [0., 100.]])
    mass_many = np.array([1., 2., 1., 3, 1., 2., 1., 3, 1., 2., 1., 3])
    dt = 1 / 30.
    t_max = 30
    radius = 5

    movie = mv.Movie_2d(fun=simulation_step, dt=dt, t_max=t_max, loc=loc_many,
                        vel=vel_many, domain=domain, mass=mass_many, radius=radius)
    movie.animate("test_2d_many_balls_new")
