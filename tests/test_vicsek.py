from Code.squirmer import Squirmer
from Code.vicsek import Vicsek_continous
import numpy as np

def test_vicsek():
    N = 20
    R = 0.25
    L = 10.0
    v0 = 1.0
    beta = 0.5
    radius = 0.1
    T = 1
    dt = 0.1
    noise = 0.1
    Es = 1

    vicsek_model = Vicsek_continous(N, R, L, v0, beta, radius, Es, T, dt, noise)
    #Runs the simulation and plot at intervals
    num_steps = 3
    for _ in range(num_steps):
        vicsek_model.loop_time()
        prct = vicsek_model.how_many_in_square()
        assert prct == 100