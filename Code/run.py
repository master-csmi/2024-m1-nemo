import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
from squirmer import Squirmer
from interactingsquirmers import InteractingSquirmers

def plot_3squirmers(R, historyb0, historybinf, historybsup, filename, dir='graphs'):
    #Plot 3 squirmers with the same initial orientation but with different beta
    plt.figure(figsize=(8, 8))
    plt.plot([-R, R], [-R, -R], 'k-', linewidth=2)
    plt.plot([-R, R], [R, R], 'k-', linewidth=2)
    plt.plot([-R, -R], [-R, R], 'k-', linewidth=2)
    plt.plot([R, R], [-R, R], 'k-', linewidth=2)

    histories = [(historyb0, 'beta=0'), (historybinf, 'beta<0'), (historybsup, 'beta>0')]
    colors = ['blue', 'green', 'red']

    for history, label in zip(histories, colors):
        squirmer1_x, squirmer1_y, squirmer1_orient = [], [], []
        squirmer2_x, squirmer2_y, squirmer2_orient = [], [], []

        for step in history[0]:
            squirmer1_x.append(step[0])
            squirmer1_y.append(step[1])
            squirmer1_orient.append(step[4])
            squirmer2_x.append(step[2])
            squirmer2_y.append(step[3])
            squirmer2_orient.append(step[5])
        
        plt.plot(squirmer1_x, squirmer1_y, label=f'Squirmer1 {history[1]}', color=label)
        plt.plot(squirmer2_x, squirmer2_y, label=f'Squirmer2 {history[1]}', color=label)

    #Plot initial orientations
    plt.quiver(squirmer2_x[0], squirmer2_y[0], np.cos(squirmer2_orient[0]), np.sin(squirmer2_orient[0]), color='black', scale=50, width=0.002)
    plt.quiver(squirmer1_x[0], squirmer1_y[0], np.cos(squirmer1_orient[0]), np.sin(squirmer1_orient[0]), color='black', scale=50, width=0.002)

    plt.axis('equal')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Positions and Orientations of Squirmers')
    plt.legend()
    plt.grid(True)

    if not os.path.exists(dir):
        os.makedirs(dir)
    save_path = os.path.join(dir, filename + '.png')
    plt.savefig(save_path)
    plt.close()

# Define parameters
v0 = 1
L = 2
R = L/2
a = 0.15
x1, y1 = -a/1.1, 0
x2, y2 = a/1.1, 0
orient1, orient2 = np.pi/2, np.pi/2
beta0 = 0
betainf = -7.5
betasup = 7.5
dt = 1e-4
lnEps_cr = np.exp(-5)
Es = 1
T = 2
dt_out = 0.05
Eo = ((3./10.)*v0/a)
ds = 2**(7./6)*a
mu = 0.01

squirmer1b0 = Squirmer(x1, y1, orient1, a, beta0, v0)
squirmer2b0 = Squirmer(x2, y2, orient2, a, beta0, v0)

squirmer1binf = Squirmer(x1, y1, orient1, a, betainf, v0)
squirmer2binf = Squirmer(x2, y2, orient2, a, betainf, v0)

squirmer1bsup = Squirmer(x1, y1, orient1, a, betasup, v0)
squirmer2bsup = Squirmer(x2, y2, orient2, a, betasup, v0)

interact_sqb0 = InteractingSquirmers(squirmer1b0, squirmer2b0, R, dt, dt_out, T, Es, ds, mu, Eo, lnEps_cr)
historyb0 = interact_sqb0.loop_time()

interact_sqbinf = InteractingSquirmers(squirmer1binf, squirmer2binf, R, dt, dt_out, T, Es, ds, mu, Eo, lnEps_cr)
historybinf = interact_sqbinf.loop_time()

interact_sqbsup = InteractingSquirmers(squirmer1bsup, squirmer2bsup, R, dt, dt_out, T, Es, ds, mu, Eo, lnEps_cr)
historybsup = interact_sqbsup.loop_time()

plot_3squirmers(R, historyb0, historybinf, historybsup, "diff_betatest")