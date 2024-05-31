import numpy as np
from squirmer import Squirmer
from interactingsquirmers import InteractingSquirmers
from plot import plot_sim_6squirmers

# Define parameters
#velocity
v0 = 1
#length of the square
L = 2
#half of the length of the square
R = L/2
#squirmers' radius
a = 0.15
#coordinates
x1, y1 = -a/1.1, 0
x2, y2 = a/1.1, 0
#orientations
orient1 = np.pi/2
orient2 = [-np.pi/4, np.pi/4, np.pi/2, -np.pi/2, 3*np.pi/4, -3*np.pi/4, np.pi, 2*np.pi]
#files labels
file = ["-pi.4", "pi.4", "pi.2", "-pi.2", "3pi.4", "-3pi.4", "pi", "2pi"]
#betas
beta0 = 0
betainf = -7.5
betasup = 7.5
#time-step
dt = 1e-4
#cut-off for -log
lnEps_cr = np.exp(-5)
#amplitude of steric interactions
Es = 1
#simulation time
T = 1.2
#periodicity of outputs
dt_out = 0.05
#amplitude of orientational interactions
Eo = ((3./10.)*v0/a)
#distance of steric interactions
ds = 2**(7./6)*a
#viscovity parameter
mu = 0.01

for orient, pi in zip(orient2, file):

    squirmer1b0 = Squirmer(x1, y1, orient1, a, beta0, v0)
    squirmer2b0 = Squirmer(x2, y2, orient, a, beta0, v0)

    squirmer1binf = Squirmer(x1, y1, orient1, a, betainf, v0)
    squirmer2binf = Squirmer(x2, y2, orient, a, betainf, v0)

    squirmer1bsup = Squirmer(x1, y1, orient1, a, betasup, v0)
    squirmer2bsup = Squirmer(x2, y2, orient, a, betasup, v0)

    interact_sqb0 = InteractingSquirmers(squirmer1b0, squirmer2b0, R, dt, dt_out, T, Es, ds, mu, Eo, lnEps_cr)
    historyb0 = interact_sqb0.loop_time()

    interact_sqbinf = InteractingSquirmers(squirmer1binf, squirmer2binf, R, dt, dt_out, T, Es, ds, mu, Eo, lnEps_cr)
    historybinf = interact_sqbinf.loop_time()

    interact_sqbsup = InteractingSquirmers(squirmer1bsup, squirmer2bsup, R, dt, dt_out, T, Es, ds, mu, Eo, lnEps_cr)
    historybsup = interact_sqbsup.loop_time()

    plot_sim_6squirmers(R, historyb0, historybinf, historybsup, f"sq2.{pi}")

