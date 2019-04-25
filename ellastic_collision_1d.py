import numpy
import matplotlib.pyplot as plt
import matplotlib.animation as anm                                       

def advection_1d(x_i, v, dt=1):
    return x_i + v * dt


def collision_1d(v1_i, v2_i, m1=1, m2=1):
    v1_f = (m1 - m2) / (m1 + m2) * v1_i + 2 * m2 / (m1 + m2) * v2_i
    v2_f = (m2 - m1) / (m1 + m2) * v2_i + 2 * m1 / (m1 + m2) * v1_i
    return v1_f, v2_f


def simulation_step(dt, x1_0, x2_0, v1_0, v2_0, domain_x):
    """ returns positions and velocity  of the particles after one step"""
    # TODO - assuming particle radius = 1
    r1 = .3
    r2 = .3
    x1 = x1_0
    x2 = x2_0
    v1 = v1_0
    v2 = v2_0

    # do collision between 2 particles
    x1_i, x2_i = x1, x2
    x1 = advection_1d(x1_i, dt, v1)
    x2 = advection_1d(x2_i, dt, v2)
    if abs(x1-x2) < (r1 + r2):
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

# global simulation parameters (should we switch to class?)
fps=30.
n_max=400
dt=1/30.
loc=[1, 5]
vel=[1, -1]
domain_x=[0,10]

# plotting
fig = plt.figure()
fig.subplots_adjust(left=0, bottom=0, right=1, top=1)

ax = fig.add_subplot(111, aspect='equal', autoscale_on=False,xlim=(domain_x), ylim=(0,1))
plot_particles, = ax.plot([], [], 'bo')

frame = plt.Rectangle([0,0], domain_x[1] - domain_x[0], 1, fc='none')
ax.add_patch(frame)

def generate_frame(i):
    global plot_particles, frame, dt, loc, vel, domain_x
    loc, vel = simulation_step(dt, loc[0], loc[1], vel[0], vel[1], domain_x)
    plot_particles.set_data(loc, .5)
    plot_particles.set_markersize(30) # TODO - set size properly
    return plot_particles, frame

# do animation and run the simulation
animation = anm.FuncAnimation(fig, generate_frame, frames=n_max, interval=10, blit=True)

# save animation
animation.save('1d_dorota_example.mp4', fps=30, extra_args=['-vcodec', 'libx264'])
