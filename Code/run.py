import numpy as np
from squirmer import Squirmer
from interactingsquirmers import InteractingSquirmers

# Define parameters
v0 = 1
L = 2
R = L/2
a = 0.15
x1, y1 = -0.7, 0.7
x2, y2 = a/1.1, 0
orient1, orient2 = 3*np.pi/4, 0
beta = 0
dt = 1e-4
lnEps_cr = np.exp(-5)
Es = 1
T = 1
dt_out = 0.05
Eo = ((3./10.)*v0/a)
ds = 2**(7./6)*a
mu = 0.01

squirmer1 = Squirmer(x1, y1, orient1, a, beta, v0)
squirmer2 = Squirmer(x2, y2, orient2, a, beta, v0)

interact_sq = InteractingSquirmers(squirmer1, squirmer2, R, dt, dt_out, T, Es, ds, mu, Eo, lnEps_cr)
interact_sq.run("bordertest", "bordertest", "dist_bord_test")
     