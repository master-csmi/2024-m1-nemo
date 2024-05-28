import numpy as np
from squirmer import Squirmer
from interactingsquirmers import InteractingSquirmers
from plot import plot_3squirmers

# Define parameters
v0 = 1
L = 2
R = L/2
a = 0.15
x1, y1 = -a/1.1, 0
x2, y2 = a/1.1, 0
orient1, orient2 = np.pi/2, np.pi
beta0 = 0
betainf = -7.5
betasup = 7.5
dt = 1e-4
lnEps_cr = np.exp(-5)
Es = 1
T = 0.7
dt_out = 0.05
Eo = ((3./10.)*v0/a)
ds = 2**(7./6)*a
mu = 0.01

squirmer1b0 = Squirmer(x1, y1, orient1, a, beta0, v0)
squirmer2b0 = Squirmer(x2, y2, orient2, a, beta0, v0)

squirmer1binf = Squirmer(x1, y1, orient1, a, betainf, v0)
squirmer2binf = Squirmer(x2, y2, orient2, a, betainf, v0)

squirmer1bsup = Squirmer(x1, y1, orient1, a, betasup, v0)
squirmer2bsup = Squirmer(x2, y2, orient2, a, betasup, v0)

interact_sqb0 = InteractingSquirmers(squirmer1b0, squirmer2b0, R, dt, dt_out, T, Es, ds, mu, Eo, lnEps_cr)
#historyb0 = interact_sqb0.loop_time()

interact_sqbinf = InteractingSquirmers(squirmer1binf, squirmer2binf, R, dt, dt_out, T, Es, ds, mu, Eo, lnEps_cr)
#historybinf = interact_sqbinf.loop_time()

interact_sqbsup = InteractingSquirmers(squirmer1bsup, squirmer2bsup, R, dt, dt_out, T, Es, ds, mu, Eo, lnEps_cr)
#historybsup = interact_sqbsup.loop_time()

#interact_sqb0.run("sq2=pi.", "sq2_pi.", "dist_sq2_pi.")
#interact_sqbinf.run("sq2=pi.", "sq2_pi.", "dist_sq2_pi.")
interact_sqbsup.run("sq2=pi.", "sq2_pi.", "dist_sq2_pi.")

#plot_3squirmers(R, historyb0, historybinf, historybsup, "sq2.-pi.4")