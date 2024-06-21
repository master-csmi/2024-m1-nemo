from squirmer import Squirmer
from interactingsquirmers import InteractingSquirmers
from plot import plot_sim_nsquirmers, plot_sim_squirmer_border, create_video_from_history

def sim_interacting_squirmers(N, xs, ys, orients, a, beta, v0, R, dt, dt_out, T, Es, ds, mu, Eo, lnEps_cr, border, sim_bord, filename, dir='graphs'):
    interact = InteractingSquirmers(N, xs, ys, orients, a, beta, v0, R, dt, dt_out, T, Es, ds, mu, Eo, lnEps_cr, border)
    history = interact.loop_time()

    plot_sim_nsquirmers(history, R, N, sim_bord, filename=filename, dir=dir)

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