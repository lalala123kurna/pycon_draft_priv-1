from pycontest import ellastic_collision_2d_many as ec                         
import numpy as np                                                             
import pytest 

def test_regression():

    # initial condition and simulation params
    vel_0 = np.array([[12., 12.], [-8, -15.], [-3, 4],     [10, 0],    [2., 12.],  [-8, -15.],
                      [3,    -4], [-10, 0],   [-12., 12.], [8, -15.],  [-3, -4],   [-10, 0]])
    loc_0 = np.array([[20., 45.], [65, 70.],  [85., 90.],  [10., 10.], [50., 45.], [15, 70.],
                      [15., 90.], [45., 10.], [75., 45.],  [40, 70.],  [45., 90.], [90., 10.]])
    domain = np.array([[0., 100.], [0., 100.]])
    mass = np.array([1., 2., 1., 3, 1., 2., 1., 3, 1., 2., 1., 3])
    dt = 1 / 30.
    t_max = 30
    radius = 5
    t=dt

    # reference values
    loc_ref = np.array([[45.08684936, 91.01680413], [55.03873045, 77.89907846],
                        [35.83978766, 82.60457335], [20.01748873, 81.17175372],
                        [ 7.04014445, 30.65715872], [80.2073723 , 50.87462212],
                        [59.82053928, 44.32842005], [ 8.04620836, 44.14156685],
                        [31.5298449 , 95.31191695], [ 9.62453512, 61.1300588 ],
                        [55.30100701,  9.8921577 ], [ 8.44782247, 86.48334328]])

    vel_ref = np.array([[ 18.23726142,  23.13062431], [  2.40950758,   6.53985751],
                        [  0.806431  ,  -9.92303699], [ -4.93232048,  -0.41164134],
                        [ -8.58411706,  26.30953988], [ -6.06871703,  13.04029963],
                        [  3.12367898, -14.5890102 ], [-12.38423597,  -1.33033359],
                        [-14.516909  ,  -4.37959354], [  0.42514493,   0.76911109],
                        [  0.29948754,   1.27290205], [  5.06255334,   6.26903091]])

    # run simulation
    loc = np.copy(loc_0)                                                       
    vel = np.copy(vel_0)                                                       
    while(t<t_max):                                                            
        loc, vel = ec.simulation_step(dt, mass, radius, loc, vel, domain)      
        t += dt                

    # from scipy docs:
    # The tolerance values are positive, typically very small numbers.
    # The relative difference (rtol * abs(b)) and the absolute difference atol
    # are added together to compare against the absolute difference 
    # between a and b.
    np.testing.assert_allclose(loc, loc_ref, atol=0, rtol=1e-9, err_msg="for locations")
    np.testing.assert_allclose(vel, vel_ref, atol=0, rtol=1e-7, err_msg="for velocities")
