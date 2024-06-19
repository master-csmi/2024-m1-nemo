import argparse
import numpy as np
from simulation import sim_interacting_squirmers, sim_vid_interact_sq, sim_sq_border

def main(simulation, N, filename):
    # Define parameters
    #velocity
    v0 = 1
    #length of the square
    L = 2
    #half of the length of the square
    R = L/2
    #squirmers' radius
    a = 0.05
    #betas
    betas = [0, -7.5, 7.5]
    beta = -7.5
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
    #viscovity parameter
    mu = 0.01
    #amplitude of orientational interactions
    Eo = [((3./10.)*v0/a), (16/10)*mu*np.pi*a**2, (-3./2.)*(v0/a)]
    #distance of steric interactions
    ds = 2**(7./6)*a

    #coordinates and orientations
    # xs, ys, orients = np.zeros(N, dtype=float), np.zeros(N, dtype=float), np.zeros(N, dtype=float)
    # for j in range(N):
    #     orientation = np.random.uniform(0, 2*np.pi)
    #     orients[j] = orientation
    # k = 0
    # while k < N:
    #     x = np.random.uniform(-(R-2*a), (R-2*a))
    #     y = np.random.uniform(-(R-2*a), (R-2*a))
    #     #Each particle must have a unique initial position
    #     if not any(np.isclose(x, xs, atol=1.1*a)) and not any(np.isclose(y, ys, atol=1.1*a)):
    #         xs[k] = x
    #         ys[k] = y
    #         k += 1
    xs = [-a, 2*a/1.5]
    ys = [0, 0]
    orients = [np.pi/2, np.pi/2]
    border = False
    border_plot = False

    # inter = InteractingSquirmers(N, xs, ys, orients, a, beta, v0, R, dt, dt_out, T, Es, ds, mu, Eo, lnEps_cr, border)
    # sim_vid_interact_sq(N, xs, ys, orients, a, beta, v0, R, dt, dt_out, T, Es, ds, mu, Eo, lnEps_cr, border, filename)
    # inter.loop_time()

    if simulation == 'video':
        sim_vid_interact_sq(N, xs, ys, orients, a, beta, v0, R, dt, dt_out, T, Es, ds, mu, Eo[0], lnEps_cr, border, filename, dir='videos')
    elif simulation == 'plot':
        sim_interacting_squirmers(N, xs, ys, orients, a, beta, v0, R, dt, dt_out, T, Es, ds, mu, Eo[0], lnEps_cr, border, filename, border_plot, dir='graphs')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run squirmer simulations.")
    parser.add_argument('simulation', choices=['video', 'plot'],
                        help="Choose which simulation to run")
    parser.add_argument('N', type=int, help="Number of squirmer")
    parser.add_argument('filename', type=str, help="Filename for saving the results")
    args = parser.parse_args()

    main(args.simulation, args.N, args.filename)
