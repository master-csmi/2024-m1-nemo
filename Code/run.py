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
    mu = 1
    #amplitude of orientational interactions
    #Eo[0] = Eoinitial, E[1] = Eobrumley, E[2] = Eolauga
    Eo = [(((3./10.)*v0/a), "Eo_init"), ((16/10)*mu*np.pi*a**2, "Eo_brumley"), ((-3./2.)*(v0/a), "Eo_lauga")]
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
    border = False
    border_plot = False

    # inter = InteractingSquirmers(N, xs, ys, orients, a, beta, v0, R, dt, dt_out, T, Es, ds, mu, Eo, lnEps_cr, border)
    # sim_vid_interact_sq(N, xs, ys, orients, a, beta, v0, R, dt, dt_out, T, Es, ds, mu, Eo, lnEps_cr, border, filename)
    # inter.loop_time()

    if simulation == 'video':
        sim_vid_interact_sq(N, xs, ys, orients, a, beta, v0, R, dt, dt_out, T, Es, ds, mu, Eo[0][0], lnEps_cr, border, filename, dir='videos')
    elif simulation == 'plot':
        sim_interacting_squirmers(N, xs, ys, orients, a, beta, v0, R, dt, dt_out, T, Es, ds, mu, Eo[0][0], lnEps_cr, border, filename, border_plot, dir='graphs')
    elif simulation == 'Eo_sim':
        betas = [(0, "beta0"), (-7.5, "betainf"), (7.5, "betasup")]
        while True:
            output_type = input("Which type of simulation? (plot/video): ").strip().lower()
            if output_type in ['plot', 'video']:
                break
            else:
                print("Invalid input. Please enter 'plot' or 'video'.")
        xseo = [-a, 2*a]
        yseo = [0, 0]
        orient1 = np.pi/2
        orient2 = [(np.pi/2, "pi_2_"), (-np.pi/2, "mpi_2_"), (3*np.pi/4, "3pi_4_"), (-3*np.pi/4, "m3pi_4_"), (np.pi, "pi_"), (2*np.pi, "2pi_"), (np.pi/4, "pi_4"), (-np.pi/4, "mpi_4")]
        for (betaeo, labelbeta) in betas:
            for (Eos, labeleo) in Eo:
                for (pi, labelpi) in orient2:
                    filenameeo = labelpi + labeleo
                    orientseo = [orient1, pi]
                    if output_type == 'plot':
                        direo = 'graphs/Eo_analysis/' + labelbeta + '/' + labelpi
                        sim_interacting_squirmers(2, xseo, yseo, orientseo, a, betaeo, v0, R, dt, dt_out, T, Es, ds, mu, Eos, lnEps_cr, border, filenameeo, border_plot, dir=direo)
                    else:
                        direo = 'videos/Eo_analysis/' + labelbeta + '/' + labelpi
                        sim_vid_interact_sq(2, xseo, yseo, orientseo, a, betaeo, v0, R, dt, dt_out, T, Es, ds, mu, Eos, lnEps_cr, border, filenameeo, dir=direo)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run squirmer simulations.")
    parser.add_argument('simulation', choices=['video', 'plot', 'Eo_sim'],
                        help="Choose which simulation to run")
    parser.add_argument('N', type=int, help="Number of squirmer")
    parser.add_argument('filename', type=str, help="Filename for saving the results")
    args = parser.parse_args()

    main(args.simulation, args.N, args.filename)
