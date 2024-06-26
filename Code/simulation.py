import numpy as np
import time
import matplotlib.pyplot as plt
from squirmer import Squirmer
from interactingsquirmers import InteractingSquirmers
from vicsek import Vicsek_continous
from plot import plot_sim_nsquirmers, plot_sim_squirmer_border, create_video_from_history, plot_time

def sim_interacting_squirmers(N, xs, ys, orients, a, beta, v0, Nx, Ny, dt, dt_out, T, Es, ds, mu, Eo, lnEps_cr, Do, no, border, sim_bord, filename, border_plot, dir='graphs'):
    #border : False or True to plot the borders
    #sim_border : False or True, if True it does the simulation with one squirmer and one border
    interact = InteractingSquirmers(N, xs, ys, orients, a, beta, v0, Nx, Ny, dt, dt_out, T, Es, ds, mu, Eo, lnEps_cr, Do, no, border)
    history = interact.loop_time()
    plot_time(interact, interact.vector_dists_min, "min_dist_" + filename, 'minimal distance', dir)
    plot_time(interact, interact.list_polar, "polar_" + filename, 'polar parameter', dir)

    plot_sim_nsquirmers(history, Nx, Ny, N, a, border_plot, sim_bord, filename=filename, dir=dir)

def sim_sq_border(x_positions, y_positions, orientations, a, betas, v0, R, dt, dt_out, T, Es, ds, mu, Eo, lnEps_cr, Do, no, filename, dir='graphs'):
    #Compute simulation for 4 squirmers close to a border
    histories = []
    for beta in betas:
        squirmers = [Squirmer(x, y, orient, a, beta, v0) for x, y, orient in zip(x_positions, y_positions, orientations)]
        interact_sqs = [InteractingSquirmers(squirmers[i], squirmers[i + 1], R, dt, dt_out, T, Es, ds, mu, Eo, lnEps_cr, Do, no) for i in range(0, len(squirmers), 2)]
        histories.append([interact_sq.loop_time() for interact_sq in interact_sqs])

    plot_sim_squirmer_border(R, histories, filename, dir=dir)

def sim_vid_interact_sq(N, xs, ys, orients, a, beta, v0, Nx, Ny, dt, dt_out, T, Es, ds, mu, Eo, lnEps_cr, Do, no, border, filename, dir='videos'):
    #Create a video of N squirmers interacting
    interact_sq = InteractingSquirmers(N, xs, ys, orients, a, beta, v0, Nx, Ny, dt, dt_out, T, Es, ds, mu, Eo, lnEps_cr, Do, no, border)
    history = interact_sq.loop_time()
    plot_time(interact_sq, interact_sq.vector_dists_min, "min_dist_" + filename, 'minimal distance', dir)
    plot_time(interact_sq, interact_sq.list_polar, "polar_" + filename, 'polar parameter', dir)

    create_video_from_history(history, Nx, Ny, N, a, filename=filename, dir=dir)

def sim_Eo_param(Eo, a, v0, dt, dt_out, T, Es, ds, mu, lnEps_cr, Do, no, border, border_plot):
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
            sim_interacting_squirmers(2, xseo, yseo, orienteo, a, betas[0][0], v0, 1, 1, dt, dt_out, T, Es, ds, mu, Eos, lnEps_cr, Do, no, border, filenameeo, border_plot, dir=direo)
        else:
            direo = 'videos/Eo_analysis/' + betas[0][1] + '/' + orient2[2][1]
            sim_vid_interact_sq(2, xseo, yseo, orienteo, a, betas[0][0], v0, 1, 1, dt, dt_out, T, Es, ds, mu, Eos, lnEps_cr, Do, no, border, filenameeo, dir=direo)

def sim_vicsek(N, R, L, v0, beta, radius, T, dt, no, nb_step):

    vicsek_model = Vicsek_continous(N, R, L, v0, beta, radius, T, dt, no)

    #Plots initial positions and save the figure
    fig, ax = plt.subplots(figsize=(8, 8))
    vicsek_model.plot(ax)
    plt.title("Initial Positions")
    plt.savefig("graphs/vicsek/vicsek_initial_positions.png")
    plt.close()
    polar = vicsek_model.polar_order_parameter()
    print(f"polar parameter = {polar}")
    #Runs the simulation and plot at intervals
    for step in range(nb_step):
        start_time = time.time()
        vicsek_model.loop_time()
        end_time = time.time()
        sim_time = end_time - start_time
        print(f"Simulation {step + 1} took {sim_time:.2f} seconds")
        prct = vicsek_model.how_many_in_square()

        polar = vicsek_model.polar_order_parameter()
        print(f"polar parameter = {polar}")

        fig, ax = plt.subplots(figsize=(8, 8))
        vicsek_model.plot(ax)
        plt.title(f"Positions at Step {step + 1}")
        plt.savefig(f"graphs/vicsek/vicsek_positions_step_{step + 1}.png")
        plt.close()