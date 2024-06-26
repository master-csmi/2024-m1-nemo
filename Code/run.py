import argparse
import numpy as np
from simulation import sim_interacting_squirmers, sim_vid_interact_sq, sim_sq_border, sim_Eo_param, sim_vicsek

def main(simulation, N, filename):
    # Define parameters
    #velocity
    v0 = 1
    #length of the square
    L = 5
    #half of the length of the square
    R = L/2
    #squirmers' radius
    a = 0.15
    #betas
    beta = 0
    #time-step
    dt = 1e-4
    #cut-off for -log
    lnEps_cr = np.exp(-5)
    #amplitude of steric interactions
    Es = 1
    #simulation time
    T = 0.5
    #periodicity of outputs
    dt_out = 0.05
    #viscovity parameter
    mu = 0.1
    #amplitude of orientational interactions
    #Eo[0] = Eoinitial, E[1] = Eobrumley, E[2] = Eolauga
    Eo = [(((3./10.)*v0/a), "Eo_init"), ((16/10)*mu*np.pi*a**2, "Eo_brumley"), ((-3./2.)*(v0/a), "Eo_lauga"),
          (((-3./10.)*v0/a), "mEo_init"), (-5, "m5"), (0.005, "m0_005"), (-2, "m2"), (-0.5, "m0_5"),
          (0.5, "0_5")]
    #distance of steric interactions
    ds = 2**(7./6)*a

    #coordinates and orientations
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
    #border to simulate chanel or box
    border = False
    #border_plot to do the simulation of the border or not
    border_plot = True

    # inter = InteractingSquirmers(N, xs, ys, orients, a, beta, v0, R, dt, dt_out, T, Es, ds, mu, Eo, lnEps_cr, border)
    # sim_vid_interact_sq(N, xs, ys, orients, a, beta, v0, R, dt, dt_out, T, Es, ds, mu, Eo, lnEps_cr, border, filename)
    # inter.loop_time()

    if simulation == 'video':
        sim_vid_interact_sq(N, xs, ys, orients, a, beta, v0, R, dt, dt_out, T, Es, ds, mu, Eo[0][0], lnEps_cr, border, filename, dir='videos')
    elif simulation == 'plot':
        xs = [-a, a/0.9]
        ys = [0, 0]
        orients = [np.pi/2, np.pi/2]
        sim_interacting_squirmers(N, xs, ys, orients, a, beta, v0, R, dt, dt_out, T, Es, ds, mu, Eo[0][0], lnEps_cr, border, False, filename, border_plot, dir='graphs')
    elif simulation == 'border':
        xs = [-0.4]
        ys = [-0.7]
        orient = [[-np.pi/6], [-np.pi/4], [-np.pi/3], [-np.pi/2]] 
        N = 1
        sim_border = True
        for i, pi in enumerate(orient):
            filename = 'sim_num_' + str(i)
            sim_interacting_squirmers(N, xs, ys, pi, a, beta, v0, R, dt, dt_out, T, Es, ds, mu, Eo[0][0], lnEps_cr, border, sim_border, filename, border_plot, dir='graphs/border')
    elif simulation == 'Eo_sim':
        sim_Eo_param(Eo, a, v0, dt, dt_out, T, Es, ds, mu, lnEps_cr, border, border_plot)
    elif simulation == 'vicsek':
        N = 20
        R = 0.25
        L = 10.0
        v0 = 1.0
        beta = 0.5
        radius = 0.1
        T = 1
        dt = 0.1
        noise = 1e-4
        nb_step = 5
        sim_vicsek(N, R, L, v0, beta, radius, T, dt, noise, nb_step)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run squirmer simulations.")

    parser.add_argument('simulation', choices=['video', 'plot', 'Eo_sim','border', 'vicsek'],
                        help="Choose which simulation to run")
    parser.add_argument('N', type=int, help="Number of squirmer")
    parser.add_argument('filename', type=str, help="Filename for saving the results")
    args = parser.parse_args()

    main(args.simulation, args.N, args.filename)
