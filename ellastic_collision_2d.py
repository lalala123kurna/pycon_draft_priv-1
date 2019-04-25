import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anm                                       
import math as mt

def advection_1d(x_i, v, dt=1):
    return x_i + v * dt


def collision_1d(v1_i, v2_i, m1=1, m2=1):
    v1_f = (m1 - m2) / (m1 + m2) * v1_i + 2 * m2 / (m1 + m2) * v2_i
    v2_f = (m2 - m1) / (m1 + m2) * v2_i + 2 * m1 / (m1 + m2) * v1_i
    return v1_f, v2_f


def collision_2d(v1, v2, dist, m1=1, m2=1):
    n = dist / np.linalg.norm(dist)
    # scalars
    v1n = np.dot(v1, n)
    v2n = np.dot(v2, n)
    v1s = v1 - v1n * n
    v2s = v2 - v2n * n

    v1n_f, v2n_f = collision_1d(v1n, v2n, m1, m2) * n
    v1_f = v1s + v1n_f
    v2_f = v2s + v2n_f
    return v1_f, v2_f


def simulation_step(dt, loc1_0, loc2_0, v1_0, v2_0, domain):
    """ returns positions and velocity  of the particles after one step"""

    # radius
    r1 = r2 = 20

    # positions  = (x, y)
    loc1 = loc1_0
    loc2 = loc2_0
    dist = loc2 - loc1

    # elocities
    v1 = v1_0
    v2 = v2_0

    # do collision between 2 particles
    for i in [0, 1]:
        loc1[i] = advection_1d(loc1[i], dt, v1[i])
        loc2[i] = advection_1d(loc2[i], dt, v2[i])
    if np.linalg.norm(dist) < (r1 + r2):
        v1, v2 = collision_2d(v1, v2, dist)

    #boundary condition
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

    def __init__(self, dt, t_max, loc1, loc2, vel1, vel2, domain, ms=30):

        self.dt = dt                                                           
        self.fps = int(mt.floor(1/dt))                                         
        self.t_max = t_max                                                     
        self.n_max = int(mt.floor(t_max / dt))                                 
        self.loc1 = loc1                                                         
        self.vel1 = vel1
        self.loc2 = loc2                                                         
        self.vel2 = vel2
        self.domain = domain                                               
        self.ms = ms

        self.fig = plt.figure()                                                
        self.fig.subplots_adjust(left=0, bottom=0, right=1, top=1)             

        self.ax = self.fig.add_subplot(111, aspect='equal', autoscale_on=False, xlim=(domain[0]), ylim=domain[1])
        self.plot_particles, = self.ax.plot([], [], 'bo')                      

        self.frame = plt.Rectangle([0,0], domain[0][1] - domain[0][0], domain[1][1] - domain[1][0], fc='none')
        self.ax.add_patch(self.frame)        

    def generate_frame(self, i):                                               
        loc_arr, vel_arr = simulation_step(self.dt, self.loc1, self.loc2, self.vel1, self.vel2, self.domain)
        self.loc1 = loc_arr[0]
        self.loc2 = loc_arr[1]
        self.vel1 = vel_arr[0]
        self.vel2 = vel_arr[1]
        self.plot_particles.set_data([self.loc1[0], self.loc2[0]], [self.loc1[1], self.loc2[1]])                              
        self.plot_particles.set_markersize(self.ms)                            
        return self.plot_particles, self.frame                                 
                                                                               
    def animate(self, name):                                                   
        animation = anm.FuncAnimation(self.fig, self.generate_frame, frames=self.n_max, interval=10, blit=True)
                                                                               
        fps = self.fps                                                         
        animation.save(name+'.mp4', fps=fps, extra_args=['-vcodec', 'libx264'])

loc1 = np.array([30, 50])
loc2 = np.array([30 + (1200) ** 0.5, 50 + 20])
v1 = [12, 12]
v2 = [-8, -15]
domain = [[0, 100], [0, 100]]

movie = Movie_2d(1/30., 10, loc1, loc2, v1, v2, domain)                             
movie.animate("test_2d")   
