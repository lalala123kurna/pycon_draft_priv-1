"""
Animation of Elastic collisions with Gravity

author: Jake Vanderplas
email: vanderplas@astro.washington.edu
website: http://jakevdp.github.com
license: BSD
Please feel free to use and modify this, but keep the above information. Thanks!
"""
import numpy as np
import math as mt
from scipy.spatial.distance import pdist, squareform

import matplotlib.pyplot as plt
import scipy.integrate as integrate
import matplotlib.animation as animation

class ParticleBox:
    """Orbits class
    
    init_state is an [N x 2] array, where N is the number of particles:
       [[x1, vx1],
        [x2, vx2],
        ...     ]

    bounds is the size of the box: [xmin, xmax]
    """
    def __init__(self,
                 init_state = [[-1, -1],
                               [1, 0.5]],
                 bounds = [-3,2],
                 size = .25):

        self.init_state = np.asarray(init_state, dtype=float)
        self.size = size
        self.state = self.init_state.copy()
        self.time_elapsed = 0
        self.bounds = bounds

    def step(self, dt):
        """step once by dt seconds"""
        self.time_elapsed += dt
        
        # update positions
        self.state[:, 0] += dt * self.state[:, 1]

        # find pairs of particles undergoing a collision
        D = squareform(pdist(self.state[:, :1]))
        ind1, ind2 = np.where(D < 2 * self.size)
        unique = (ind1 < ind2)
        ind1 = ind1[unique]
        ind2 = ind2[unique]

        # update velocities of colliding pairs
        for i1, i2 in zip(ind1, ind2):
            # mass
            m1 = 0.05
            m2 = 0.05

            # location vector
            r1 = self.state[i1, 0]
            r2 = self.state[i2, 0]

            # velocity vector
            v1 = self.state[i1, 1]
            v2 = self.state[i2, 1]

            # TODO - why this is not working?
            # relative location & velocity vectors
            #r_rel = r1 - r2
            #v_rel = v1 - v2
            # momentum vector of the center of mass
            #v_cm = (m1 * v1 + m2 * v2) / (m1 + m2)
            # collisions of spheres reflect v_rel over r_rel
            #rr_rel = np.dot(r_rel, r_rel)
            #vr_rel = np.dot(v_rel, r_rel)
            #v_rel = 2 * r_rel * vr_rel / rr_rel - v_rel
            # assign new velocities
            #self.state[i1, 1] = v_cm + v_rel * m2 / (m1 + m2)
            #self.state[i2, 1] = v_cm - v_rel * m1 / (m1 + m2)
    
            self.state[i1, 1] = 2 * m2 / (m1 + m2) * v2
            self.state[i2, 1] = 2 * m1 / (m1 + m2) * v1

        # check for crossing boundary
        crossed_x1 = (self.state[:, 0] < self.bounds[0] + self.size)
        crossed_x2 = (self.state[:, 0] > self.bounds[1] - self.size)

        self.state[crossed_x1, 0] = self.bounds[0] + self.size
        self.state[crossed_x2, 0] = self.bounds[1] - self.size

        self.state[crossed_x1 | crossed_x2, 1] *= -1


#------------------------------------------------------------
# init state
box = ParticleBox()
dt = 1. / 30 # 30fps
#dt = 10.#TODO
n_max = 600

#------------------------------------------------------------
# set up figure and animation
fig = plt.figure()
fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
ax = fig.add_subplot(111, aspect='equal', autoscale_on=False,
                     xlim=(box.bounds[0], box.bounds[1]), ylim=(0, 1))

# particles holds the locations of the particles
particles, = ax.plot([], [], 'bo', ms=1)

# rect is the box edge
rect = plt.Rectangle([box.bounds[0], 0], #TODO 
                     box.bounds[1] - box.bounds[0],
                     1,
                     ec='none', lw=2, fc='none')
ax.add_patch(rect)

def init():
    """initialize animation"""
    global box, rect, tmp_it
    particles.set_data([], [])
    rect.set_edgecolor('none')
    tmp_it = 0
    return particles, rect

def animate(i):
    """perform animation step"""
    global box, rect, dt, ax, fig, tmp_it
    box.step(dt)
    #print(tmp_it)
    tmp_it += 1
    #print(" ")
    #print("fig dpi      ", fig.dpi)
    #print("box size     ", box.size)
    #print("fig width    ", fig.get_figwidth())
    #print("ax xbound    ", ax.get_xbound())
    #print("diff(xbound) ", np.diff(ax.get_xbound())[0])
    # TODO - how to set the size better?
    ms = int(fig.dpi * 2 * box.size
             * fig.get_figwidth() / np.diff(ax.get_xbound())[0]
            )
    #print("ms", ms)
    ms -= 20

    # update pieces of the animation
    rect.set_edgecolor('k')
    particles.set_data(box.state[:, 0], .5)
    particles.set_markersize(ms)

    #print(box.state)

    return particles, rect

ani = animation.FuncAnimation(fig, animate, frames=n_max,
                              interval=10, blit=True, init_func=init)


# save the animation as an mp4.  This requires ffmpeg or mencoder to be
# installed.  The extra_args ensure that the x264 codec is used, so that
# the video can be embedded in html5.  You may need to adjust this for
# your system: for more information, see
# http://matplotlib.sourceforge.net/api/animation_api.html
ani.save('1d_particle_box.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

#plt.show()
