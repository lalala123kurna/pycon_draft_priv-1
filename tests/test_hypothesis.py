from pycontest import ellastic_collision_2d_many as ec

import numpy as np

import pytest

# https://hypothesis.readthedocs.io/en/latest/
from hypothesis import given, strategies as st

#TODO - move energy and momentum to separate file
def E_kin(vel, mass):
    """ calculate the kinematic energy of all particles """
    vel = np.array(vel)
    mass = np.array(mass)
    return 0.5 * np.sum(mass * vel**2)

# hypothesis simple example
@given(mass1  = st.floats(min_value=.1, max_value=1e3),
       mass2  = st.floats(min_value=.1, max_value=1e3))
def test_energy_hypothesis(mass1, mass2):

    # initial condition and simulation parameters
    domain = ([0, 20], [0, 20])
    dt = 0.5
    t_max = 6
    t = 0
    loc_0 = np.array([[3, 4],[15, 2]])
    vel_0 = np.array([[1, 0.5], [-1, -.25]])
    radius = 1

    # mass randomly chosen by hypothesis
    mass = [mass1, mass2]

    # run the simulation
    loc = np.copy(loc_0)
    vel = np.copy(vel_0)
    while(t<t_max):
        loc, vel = ec.simulation_step(dt, mass, radius, loc, vel, domain)
        t += dt

    E_ini = E_kin(vel_0, mass)
    E_end = E_kin(vel, mass)

    print("testing for mass = [", mass1, ",", mass2,"]. E_end = ", E_end)

    assert E_ini == E_end
