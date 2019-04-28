import matplotlib.pyplot as plt
import matplotlib.animation as anm
import math as mt
import pdb

class Movie_2d:
    """ TODO """

    def __init__(self, fun, dt, t_max, loc, vel, domain, mass, radius=10):
        self.dt = dt
        self.fps = int(mt.floor(1 / dt))
        self.n_max = int(mt.floor(t_max / dt))
        self.loc = loc
        self.vel = vel
        self.domain = domain
        self.mass = mass
        self.radius = radius
        self.fun = fun

        self.fig = plt.figure()
        self.fig.subplots_adjust(left=0.1, bottom=0.1, right=.9, top=.9)

        self.ax = self.fig.add_subplot(111, aspect='equal', autoscale_on=False, xlim=(domain[0]), ylim=domain[1])
        for el in range(loc.shape[0]):
            self.circ = plt.Circle(self.loc[el, :], self.radius, color='r', fill=False, linewidth=2)
            self.ax.add_patch(self.circ)

        self.frame = plt.Rectangle([0, 0], domain[0][1] - domain[0][0], domain[1][1] - domain[1][0], fc='none')
        self.ax.add_patch(self.frame)

    def generate_frame(self, i):
       # pdb.set_trace()
        loc, vel = self.fun(self.dt, self.mass, self.radius, self.loc, self.vel, self.domain)
        #pdb.set_trace()
        for el in range(loc.shape[0]):
            self.circ.set_center(loc[el, :])
        return self.circ, self.frame

    def animate(self, name):
        animation = anm.FuncAnimation(self.fig, self.generate_frame, frames=self.n_max, interval=10, blit=True)

        fps = self.fps
        animation.save(name + '.mp4', fps=fps, extra_args=['-vcodec', 'libx264'])
