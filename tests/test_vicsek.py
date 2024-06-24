import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'Code')))

from vicsek import Vicsek_continous
import numpy as np

def test_vicsek():
    N = 40
    R = 0.25
    L = 10.0
    v0 = 1.0
    beta = 0.5
    radius = 0.1
    T = 2.5
    dt = 0.1
    noise = 0.1

    vicsek_model = Vicsek_continous(N, R, L, v0, beta, radius, T, dt, noise)
    #Runs the simulation and plot at intervals
    num_steps = 3
    for _ in range(num_steps):
        vicsek_model.loop_time()
        prct = vicsek_model.how_many_in_square()
        assert prct == 100