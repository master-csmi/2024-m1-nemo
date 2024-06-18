import argparse
import numpy as np
from squirmer import Squirmer
from interactingsquirmers import InteractingSquirmers
from plot import plot_squirmers_positions
from simulation import sim_intercating_squirmers, sim_vid_interact_sq, sim_sq_border

def main(simulation, filename):
    # Define parameters
    #velocity
    v0 = 1
    #length of the square
    L = 2
    #half of the length of the square
    R = L/2
    #squirmers' radius
    a = 0.1
    #coordinates
    x1, y1 = -2*a, 0.6
    x2, y2 = a/1.1, 0
    #orientations
    orient1 = np.pi
    orient2 = [-np.pi/4, np.pi/4, np.pi/2, -np.pi/2, 3*np.pi/4, -3*np.pi/4, np.pi, 2*np.pi]
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

    N = 4
    xs, ys, orients = np.zeros(N, dtype=float), np.zeros(N, dtype=float), np.zeros(N, dtype=float)
    for j in range(N):
        orientation = np.random.uniform(0, 2*np.pi)
        orients[j] = orientation
    k = 0
    while k < N:
        x = np.random.uniform(-(R-2*a), (R-2*a))
        y = np.random.uniform(-(R-2*a), (R-2*a))
        #Each particle must have a unique initial position
        if not any(np.isclose(x, xs, atol=1.1*a)) and not any(np.isclose(y, ys, atol=1.1*a)):
            xs[k] = x
            ys[k] = y
            k += 1
    border = False

    inter = InteractingSquirmers(N, xs, ys, orients, a, beta0, v0, R, dt, dt_out, T, Es, ds, mu, Eo, lnEps_cr, border)
    sim_vid_interact_sq(N, xs, ys, orients, a, beta0, v0, R, dt, dt_out, T, Es, ds, mu, Eo, lnEps_cr, filename, border)
    inter.loop_time()

    # squirmer1 = Squirmer(x1, y1, orient1, a, beta0, v0)
    # squirmer2 = Squirmer(x2, y2, orient2[6], a, beta0, v0)
    # intera = InteractingSquirmers(squirmer1, squirmer2, R, dt, dt_out, T, Es, ds, mu, Eo, lnEps_cr)
    # history = intera.loop_time()

    # plot_squirmers_positions(R, history, 'check')
    # if simulation == 'interact_sq':
    #     sim_intercating_squirmers(x1, y1, x2, y2, orient1, orient2, a, beta0, betainf, betasup, v0, R, dt, dt_out, T, Es, ds, mu, Eo, lnEps_cr, file, dir='graphs')
    # elif simulation == 'vid_interact_sq':
    #     sim_vid_interact_sq(x1, y1, x2, y2, orient1, orient2[0], a, betas[0], v0, R, dt, dt_out, T, Es, ds, mu, Eo, lnEps_cr, filename, dir='videos')
    # elif simulation == 'sq_border':
    #     sim_sq_border(xs, ys, orients, a, betas, v0, R, dt, dt_out, T, Es, ds, mu, Eo, lnEps_cr, filename, dir='graphs')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run squirmer simulations.")
    parser.add_argument('simulation', choices=['interact_sq', 'vid_interact_sq', 'sq_border'],
                        help="Choose which simulation to run")
    parser.add_argument('filename', type=str, help="Filename for saving the results")
    args = parser.parse_args()

    main(args.simulation, args.filename)
