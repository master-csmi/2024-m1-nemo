import numpy as np
from squirmer import Squirmer
from interactingsquirmers import InteractingSquirmers
from plot import plot_sim_6squirmers, plot_sim_squirmer_border

# Define parameters
v0 = 1
L = 2
R = L/2
a = 0.15
x1, y1 = -0.6, -0.2
x2, y2 = 0.6, -0.2
x3, y3 = 0.2, 0.6
x4, y4 = -0.2, -0.6
orient1, orient2, orient3, orient4 = 3*np.pi/4, np.pi/4, 3*np.pi/4, -np.pi/4
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
squirmer3b0 = Squirmer(x3, y3, orient3, a, beta0, v0)
squirmer4b0 = Squirmer(x4, y4, orient4, a, beta0, v0)

squirmer1binf = Squirmer(x1, y1, orient1, a, betainf, v0)
squirmer2binf = Squirmer(x2, y2, orient2, a, betainf, v0)
squirmer3binf = Squirmer(x3, y3, orient3, a, betainf, v0)
squirmer4binf = Squirmer(x4, y4, orient4, a, betainf, v0)

squirmer1bsup = Squirmer(x1, y1, orient1, a, betasup, v0)
squirmer2bsup = Squirmer(x2, y2, orient2, a, betasup, v0)
squirmer3bsup = Squirmer(x3, y3, orient3, a, betasup, v0)
squirmer4bsup = Squirmer(x4, y4, orient4, a, betasup, v0)

interact_sqb01 = InteractingSquirmers(squirmer1b0, squirmer2b0, R, dt, dt_out, T, Es, ds, mu, Eo, lnEps_cr)
historyb01 = interact_sqb01.loop_time()
interact_sqb02 = InteractingSquirmers(squirmer3b0, squirmer4b0, R, dt, dt_out, T, Es, ds, mu, Eo, lnEps_cr)
historyb02 = interact_sqb02.loop_time()

interact_sqbinf1 = InteractingSquirmers(squirmer1binf, squirmer2binf, R, dt, dt_out, T, Es, ds, mu, Eo, lnEps_cr)
historybinf1 = interact_sqbinf1.loop_time()
interact_sqbinf2 = InteractingSquirmers(squirmer3binf, squirmer4binf, R, dt, dt_out, T, Es, ds, mu, Eo, lnEps_cr)
historybinf2 = interact_sqbinf2.loop_time()

interact_sqbsup1 = InteractingSquirmers(squirmer1bsup, squirmer2bsup, R, dt, dt_out, T, Es, ds, mu, Eo, lnEps_cr)
historybsup1 = interact_sqbsup1.loop_time()
interact_sqbsup2 = InteractingSquirmers(squirmer3bsup, squirmer4bsup, R, dt, dt_out, T, Es, ds, mu, Eo, lnEps_cr)
historybsup2 = interact_sqbsup2.loop_time()

histories = [[historyb01, historyb02], [historybinf1, historybinf2], [historybsup1, historybsup2]]

plot_sim_squirmer_border(R, histories, "sim_sq_border1")

#plot_sim_6squirmers(R, historyb0, historybinf, historybsup, "sq2.2pi")