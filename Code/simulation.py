import numpy as np
from squirmer import Squirmer
from interactingsquirmers import InteractingSquirmers
from plot import plot_sim_6squirmers, plot_sim_squirmer_border, create_video_from_history

def sim_intercating_squirmers(x1, y1, x2, y2, orient1, orient2, a, beta0, betainf, betasup, v0,  R, dt, dt_out, T, Es, ds, mu, Eo, lnEps_cr, file, dir='graphs'):
    #Compute simulation for two squirmers for each orientation in orient2
    #orient2 and file need to have the same length
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

        plot_sim_6squirmers(R, historyb0, historybinf, historybsup, f"sq2.{pi}", dir=dir)

def sim_squirmer_border(x_positions, y_positions, orientations, a, betas, v0, R, dt, dt_out, T, Es, ds, mu, Eo, lnEps_cr, filename, dir='graphs'):
    #Compute simulation for 4 squirmers close to a border
    histories = []
    for beta in betas:
        squirmers = [Squirmer(x, y, orient, a, beta, v0) for x, y, orient in zip(x_positions, y_positions, orientations)]
        interact_sqs = [InteractingSquirmers(squirmers[i], squirmers[i + 1], R, dt, dt_out, T, Es, ds, mu, Eo, lnEps_cr) for i in range(0, len(squirmers), 2)]
        histories.append([interact_sq.loop_time() for interact_sq in interact_sqs])

    plot_sim_squirmer_border(R, histories, filename, dir=dir)

def sim_vid_interact_squir(x1, y1, x2, y2, orient1, orient2, a, beta, v0, R, dt, dt_out, T, Es, ds, mu, Eo, lnEps_cr, filename, dir='videos'):
    #Create a video of two squirmers interacting
    squirmer1 = Squirmer(x1, y1, orient1, a, beta, v0)
    squirmer2 = Squirmer(x2, y2, orient2, a, beta, v0)

    interact_sq = InteractingSquirmers(squirmer1, squirmer2, R, dt, dt_out, T, Es, ds, mu, Eo, lnEps_cr)
    history = interact_sq.loop_time()

    create_video_from_history(history, R=R, filename=filename, dir=dir)