import numpy as np
from squirmer import Squirmer
from interactingsquirmers import InteractingSquirmers



v0 = 1
L = 2
R = L/2
x1, y1 = -0.2, 0.01
x2, y2 = 0.2, 0
orient1, orient2 = 2*np.pi, np.pi
a = 0.05
beta = -7.5
dt = 1e-4
lnEps_cr = 5
Es = 1
T = 1
dt_out = 0.05
Eo = (3/10)*v0/a
ds = 2**(7/6)*a
D = 0
squirmer1 = Squirmer(x1,y1,orient1,a,beta,v0)
squirmer2 = Squirmer(x2,y2,orient2,a,beta,v0)

interact_sq = InteractingSquirmers(squirmer1, squirmer2, R, dt, dt_out, T, Es, ds, Eo, lnEps_cr)
interact_sq.run("test")