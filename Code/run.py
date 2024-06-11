import numpy as np
from simulation import sim_intercating_squirmers, sim_vid_interact_squir, sim_squirmer_border

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
#border simulation
xs = [-0.6, 0.6, 0.2, -0.2]
ys = [-0.2, -0.2, 0.6, -0.6]
orients = [3*np.pi/4, np.pi/4, 3*np.pi/4, -np.pi/4]
#files labels
file = ["-pi.4", "pi.4", "pi.2", "-pi.2", "3pi.4", "-3pi.4", "pi", "2pi"]
#betas
betas = [0, -7.5, 7.5]
beta0, betainf, betasup = 0, -7.5, 7.5
#time-step
dt = 1e-4
#cut-off for -log
lnEps_cr = np.exp(-5)
#amplitude of steric interactions
Es = 1
#simulation time
T = 0.7
#periodicity of outputs
dt_out = 0.05
#amplitude of orientational interactions
Eo = ((3./10.)*v0/a)
#distance of steric interactions
ds = 2**(7./6)*a
#viscovity parameter
mu = 0.01

sim_squirmer_border(xs, ys, orients, a, betas, v0,  R, dt, dt_out, T, Es, ds, mu, Eo, lnEps_cr, 'bordersimtest', dir='graphs')

#sim_intercating_squirmers(x1, y1, x2, y2, orient1, orient2, a, beta0, betainf, betasup, v0,  R, dt, dt_out, T, Es, ds, mu, Eo, lnEps_cr, file, dir='graphs')

#sim_vid_interact_squir(x1, y1, x2, y2, orient1, orients[0], a, betas[0], v0, R, dt, dt_out, T, Es, ds, mu, Eo, lnEps_cr, 'sim_interact_sq', dir='videos', gif=True)