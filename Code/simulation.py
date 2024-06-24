import numpy as np
from squirmer import Squirmer
from interactingsquirmers import InteractingSquirmers
from plot import plot_sim_nsquirmers, plot_sim_squirmer_border, create_video_from_history

def sim_interacting_squirmers(N, xs, ys, orients, a, beta, v0, R, dt, dt_out, T, Es, ds, mu, Eo, lnEps_cr, border, sim_bord, filename, border_plot, dir='graphs'):
    #border : False or True to plot the borders
    #sim_border : False or True, if True it does the simulation with one squirmer and one border
    interact = InteractingSquirmers(N, xs, ys, orients, a, beta, v0, R, dt, dt_out, T, Es, ds, mu, Eo, lnEps_cr, border)
    history = interact.loop_time()

    plot_sim_nsquirmers(history, R, N, a, border_plot, sim_bord, filename=filename, dir=dir)

def sim_sq_border(x_positions, y_positions, orientations, a, betas, v0, R, dt, dt_out, T, Es, ds, mu, Eo, lnEps_cr, filename, dir='graphs'):
    #Compute simulation for 4 squirmers close to a border
    histories = []
    for beta in betas:
        squirmers = [Squirmer(x, y, orient, a, beta, v0) for x, y, orient in zip(x_positions, y_positions, orientations)]
        interact_sqs = [InteractingSquirmers(squirmers[i], squirmers[i + 1], R, dt, dt_out, T, Es, ds, mu, Eo, lnEps_cr) for i in range(0, len(squirmers), 2)]
        histories.append([interact_sq.loop_time() for interact_sq in interact_sqs])

    plot_sim_squirmer_border(R, histories, filename, dir=dir)

def sim_vid_interact_sq(N, xs, ys, orients, a, beta, v0, R, dt, dt_out, T, Es, ds, mu, Eo, lnEps_cr, border, filename, dir='videos'):
    #Create a video of N squirmers interacting
    interact_sq = InteractingSquirmers(N, xs, ys, orients, a, beta, v0, R, dt, dt_out, T, Es, ds, mu, Eo, lnEps_cr, border)
    history = interact_sq.loop_time()

    create_video_from_history(history, R=R, N=N, filename=filename, dir=dir)

def sim_Eo_param(Eo, a, v0, dt, dt_out, T, Es, ds, mu, lnEps_cr, border, border_plot):
    betas = [(0, "beta0"), (-7.5, "betainf"), (7.5, "betasup")]
    while True:
        output_type = input("Which type of simulation? (plot/video): ").strip().lower()
        if output_type in ['plot', 'video']:
            break
        else:
            print("Invalid input. Please enter 'plot' or 'video'.")
    xseo = [-a, 2*a/1.5]
    yseo = [0, 0]
    orient1 = np.pi/2
    orient2 = [(np.pi/2, "pi_2_"), (-np.pi/2, "mpi_2_"), (3*np.pi/4, "3pi_4_"), (-3*np.pi/4, "m3pi_4_"), (np.pi, "pi_"), (2*np.pi, "2pi_"), (np.pi/4, "pi_4_"), (-np.pi/4, "mpi_4_")]
    # for (betaeo, labelbeta) in betas:
    for (Eos, labeleo) in Eo:
        # for (pi, labelpi) in orient2:
        #     filenameeo = labelpi + labeleo
        #     orientseo = [orient1, pi]
        filenameeo = orient2[2][1] + labeleo
        orienteo = [orient1, orient2[2][0]]
        if output_type == 'plot':
            direo = 'graphs/Eo_analysis/' + betas[0][1] + '/' + orient2[2][1]
            sim_interacting_squirmers(2, xseo, yseo, orienteo, a, betas[0][0], v0, 0.8, dt, dt_out, T, Es, ds, mu, Eos, lnEps_cr, border, filenameeo, border_plot, dir=direo)
        else:
            direo = 'videos/Eo_analysis/' + betas[0][1] + '/' + orient2[2][1]
            sim_vid_interact_sq(2, xseo, yseo, orienteo, a, betas[0][0], v0, 0.8, dt, dt_out, T, Es, ds, mu, Eos, lnEps_cr, border, filenameeo, dir=direo)