import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anm                                       

def advection_1d(x_i, v, dt=1):
    return x_i + v * dt


def collision_1d(v1_i, v2_i, m1=1, m2=1):
    v1_f = (m1 - m2) / (m1 + m2) * v1_i + 2 * m2 / (m1 + m2) * v2_i
    v2_f = (m2 - m1) / (m1 + m2) * v2_i + 2 * m1 / (m1 + m2) * v1_i
    return v1_f, v2_f


def collision_2d(v1, v2, r12, m1, m2):
    n = r12 / np.linalg.norm(r12)
    # scalars
    v1n = np.dot(v1, n)
    v2n = np.dot(v2, n)
    v1s = v1 - v1n * n
    v2s = v2 - v2n * n

    v1n_f, v2n_f = collision_1d(v1n, v2n, m1, m2) * n
    v1_f = v1s + v1n_f
    v2_f = v2s + v2n_f
    return v1_f, v2_f


def simulation_step(dt, r1_0, r2_0, v1_0, v2_0, domain):
    """ returns positions and velocity  of the particles after one step"""

    # masses
    m1 = 1
    m2 = 2
    # radius
    rad1 = rad2 = 20

    # positions  = (x, y)
    r1 = r1_0
    r2 = r2_0
    r12 = r2 - r1

    v1 = v1_0
    v2 = v2_0

    # do collision between 2 particles
    r1_i = r1
    r2_i = r2
    for i in [0, 1]:
        r1[i] = advection_1d(r1_i[i], dt, v1[i])
        r2[i] = advection_1d(r2_i[i], dt, v2[i])

    if np.linalg.norm(r12) < (rad1 + rad2):
        v1, v2 = collision_2d(v1, v2, r12, m1, m2)

    #boundary condition
    for i in [0, 1]:
        if r1[i] - rad1 < domain[i][0]:
            v1[i] *= -1
        if r1[i] + rad1 > domain[i][1]:
            v1[i] *= -1
        if r2[i] - rad2 < domain[i][0]:
            v2[i] *= -1
        if r2[i] + rad2 > domain[i][1]:
            v2[i] *= -1

    return [r1, r2], [v1, v2]

# global simulation parameters (should we switch to class?)
fps=30.
n_max=400
dt=1/30.
r1 = np.array([30, 50])
r2 = np.array([30 + (1200) ** 0.5, 50 + 20])
v1 = [12, 12]
v2 = [-8, -15]
domain = [[0, 100], [0, 100]]

# plotting
fig = plt.figure()
fig.subplots_adjust(left=0, bottom=0, right=1, top=1)

ax = fig.add_subplot(111, aspect='equal', autoscale_on=False,
                     xlim=(domain[0]), ylim=(domain[1]))
plot_particles, = ax.plot([], [], 'bo')

frame = plt.Rectangle([0,0], domain[0][1] - domain[0][0],
                      domain[1][1] - domain[1][0], fc='none')
ax.add_patch(frame)

def generate_frame(i):
    global plot_particles, frame, dt, loc, vel, domain
    loc, vel = simulation_step(dt, r1, r2, v1, v2, domain)
    # TODO: zle przekazuje loc do set_dat
    x_loc = [loc[0][0], loc[1][0]]
    y_loc = [loc[0][1], loc[1][1]]
    plot_particles.set_data(x_loc, y_loc, .5)
    plot_particles.set_markersize(30) # TODO - set size properly
    return plot_particles, frame

# do animation and run the simulation
animation = anm.FuncAnimation(fig, generate_frame, frames=n_max, interval=10, blit=True)

# save animation
animation.save('2d_dorota_example.mp4', fps=30, extra_args=['-vcodec', 'libx264'])
