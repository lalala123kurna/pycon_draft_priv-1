from pycontest import elastic_collisions as ec
from pycontest import simulation_2d_many as sim2d
from pycontest.movies import Movie_2d
from pycontest.utils import momentum, E_kin

import numpy as np
import pytest


# simple pytest examples
def test_collision_1d_1():
    v1_f, v2_f = ec.collision_1d(v1_i=1, v2_i=-2)
    assert v1_f == -2
    assert v2_f == 1

def test_collision_1d_2():
    v1_f, v2_f = ec.collision_1d(v1_i=1, v2_i=-2, m1=2, m2=2)
    assert v1_f == -2
    assert v2_f == 1

def test_collision_1d_3():
    v1_f, v2_f = ec.collision_1d(v1_i=1, v2_i=-2, m1=1, m2=1e6)
    assert v2_f == pytest.approx(-2, rel=1e-3)


# parametrize to show how to pass test arguments
@pytest.mark.parametrize("dt", [0.5, 0.1, 0.01])
def test_simulation_1d(dt):
    # initial condition and simulation parameters
    domain = ([-2, 12], [0, 3])
    dt = dt
    t_max = 6
    loc_0 = np.array([[0, 1.5],[10, 1.5]])
    vel_0 = np.array([[1, 0], [-1, 0]])
    radius = 1
    mass = [1, 1]

    loc, vel = sim2d.simulation(t_max, dt, mass, radius, loc_0, vel_0, domain)

    # test location and velocities after colision
    assert loc[0][0] < 5
    assert loc[1][0] > 5
    assert (loc[0][1], loc[1][1]) == (loc_0[0][1], loc_0[1][1]) 

    assert vel[0][0] == -1
    assert vel[1][0] == 1
    assert (vel[0][1], vel[1][1]) == (vel_0[0][1], vel_0[1][1]) 

    # create movie
    movie = Movie_2d(sim2d.simulation_step, dt, t_max - dt, loc, vel, domain, mass, radius)
    movie.animate("pytest_movie_1d_dt_"+str(dt))


# module to show how to run simulation once and test many things
@pytest.fixture(scope="module")
def data(request):

    print("\n create data")
    # initial condition and simulation parameters
    domain = ([-2, 12], [0, 3])
    dt = 0.5
    t_max = 6
    loc_0 = np.array([[0, 1.5],[10, 1.5]])
    vel_0 = np.array([[1, 0], [-1, 0]])
    radius = 1
    mass = [1, 1]

    loc, vel = sim2d.simulation(t_max, dt, mass, radius, loc_0, vel_0, domain)

    my_data = {}
    my_data["loc_0"] = loc_0
    my_data["vel_0"] = vel_0
    my_data["loc"] = loc
    my_data["vel"] = vel
    my_data["mass"] = mass

    def data_cleanup():
        print("\n removing data")
        my_data.clear()

    request.addfinalizer(data_cleanup)                                  
    return my_data

# TODO dodac dluzszy czas
def test_energy(data):

    print("\n test energy")

    E_ini = E_kin(data["vel_0"], data["mass"])
    E_end = E_kin(data["vel"], data["mass"])

    assert E_ini == E_end

def test_momentum(data):

    print("\n test momentum")

    p_ini = momentum(data["vel_0"], data["mass"])
    p_end = momentum(data["vel"], data["mass"])

    assert np.all(p_ini == p_end)


def test_simulation_1d_fail():
    # initial condition and simulation parameters
    domain = ([-2, 12], [0, 3])
    dt = 1.
    t_max = 12
    loc_0 = np.array([[0, 1.5],[10, 1.5]])
    vel_0 = np.array([[1, 0], [-1, 0]])
    radius = 1
    mass = [1, 1]

    with pytest.raises(Exception) as excinfo:
        loc, vel = sim2d.simulation(t_max, dt, mass, radius, loc_0, vel_0, domain)
    assert "two balls are exactly in the same place" in str(excinfo.value)


#TODO - bring those up to date
# but also leave one to show mark.skip
@pytest.mark.skip(reason=" this will be infinite loop, del_x should be better defined")
def test_simulation_collision_2():
    x1_f, x2_f, time_f, v1_f, v2_f = \
        ec_1d.simulation_collision(x1_0=0, x2_0=10, v1=3, v2=-1)


@pytest.mark.skip(reason=" this will be infinite loop, we're not checking if they can meet")
def test_simulation_collision_3():
    x1_f, x2_f, time_f, v1_f, v2_f = \
        ec_1d.simulation_collision(x1_0=0, x2_0=10, v1=1, v2=1)
