import matplotlib
import numpy as np
import os
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

def plot_squirmers_positions(R, history, filename, dir='graphs'):
    #Plot the position of each squirmers
    plt.figure(figsize=(8, 8))
    plt.plot([-R, R], [-R, -R], 'k-', linewidth=2)
    plt.plot([-R, R], [R, R], 'k-', linewidth=2)
    plt.plot([-R, -R], [-R, R], 'k-', linewidth=2)
    plt.plot([R, R], [-R, R], 'k-', linewidth=2)

    squirmer1_x = []
    squirmer1_y = []
    squirmer1_orient = []
    squirmer2_x = []
    squirmer2_y = []
    squirmer2_orient = []
    time = []

    for step in history:
        squirmer1_x.append(step[0])
        squirmer1_y.append(step[1])
        squirmer1_orient.append(step[4])
        squirmer2_x.append(step[2])
        squirmer2_y.append(step[3])
        squirmer2_orient.append(step[5])
        time.append(step[-1])
    
    #Squirmers
    plt.scatter(squirmer1_x, squirmer1_y, color='blue', s=10, label = 'Squirmer 1')
    plt.scatter(squirmer2_x, squirmer2_y, color='red', s=10, label= 'Squirmer 2')
    
    for i in range(len(squirmer2_orient)):
        #Orientation
        plt.quiver(squirmer2_x[i], squirmer2_y[i], np.cos(squirmer2_orient[i]), np.sin(squirmer2_orient[i]), color='red', scale=20, width=0.002)
        plt.quiver(squirmer1_x[i], squirmer1_y[i], np.cos(squirmer1_orient[i]), np.sin(squirmer1_orient[i]), color='blue', scale=20, width=0.002)

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

def plot_dist_sq(history, filename, dir='graphs'):
    #Plot the distance between two squirmer over time
    plt.figure(figsize=(8, 6))
    dist_list = []
    time_list = []
    for step in history:
        dist_list.append(step[10])
        time_list.append(step[11])

    plt.plot(time_list, dist_list, label="Distance")
    
    plt.xlabel('Time')
    plt.ylabel('Distance between squirmers')
    plt.title('Distance between squirmers over time')
    plt.grid(True)

    if not os.path.exists(dir):
        os.makedirs(dir)
    save_path = os.path.join(dir, filename + '.png')
    plt.savefig(save_path)

def plot_sim_6squirmers(R, historyb0, historybinf, historybsup, filename, dir='graphs'):
    #Plot 2 squirmers with the same initial orientation but with different beta
    plt.figure(figsize=(8, 8))
    plt.plot([-R, R], [-R, -R], 'k-', linewidth=2)
    plt.plot([-R, R], [R, R], 'k-', linewidth=2)
    plt.plot([-R, -R], [-R, R], 'k-', linewidth=2)
    plt.plot([R, R], [-R, R], 'k-', linewidth=2)

    histories = [(historyb0, 'beta=0'), (historybinf, 'beta<0'), (historybsup, 'beta>0')]
    colors = [('blue','cyan'), ('orange','gold'), ('green','lime')]

    for history, (color1,color2) in zip(histories, colors):
        squirmer1_x, squirmer1_y, squirmer1_orient = [], [], []
        squirmer2_x, squirmer2_y, squirmer2_orient = [], [], []

        for step in history[0]:
            squirmer1_x.append(step[0])
            squirmer1_y.append(step[1])
            squirmer1_orient.append(step[4])
            squirmer2_x.append(step[2])
            squirmer2_y.append(step[3])
            squirmer2_orient.append(step[5])
        
        plt.plot(squirmer1_x, squirmer1_y, label=f'Squirmer1 {history[1]}', color=color1)
        plt.plot(squirmer2_x, squirmer2_y, label=f'Squirmer2 {history[1]}', color=color2)

        last_orient1 = squirmer1_orient[0]
        last_orient2 = squirmer2_orient[0]

        for i in range(len(squirmer1_orient)):
            if squirmer1_orient[i] != last_orient1:
                plt.quiver(squirmer1_x[i], squirmer1_y[i], np.cos(squirmer1_orient[i]), np.sin(squirmer1_orient[i]), color=color1, scale=25, width=0.005)
                last_orient1 = squirmer1_orient[i]
            if squirmer2_orient[i] != last_orient2:
                plt.quiver(squirmer2_x[i], squirmer2_y[i], np.cos(squirmer2_orient[i]), np.sin(squirmer2_orient[i]), color=color2, scale=25, width=0.005)
                last_orient2 = squirmer2_orient[i]

    #Plot initial orientations
    plt.quiver(squirmer2_x[0], squirmer2_y[0], np.cos(squirmer2_orient[0]), np.sin(squirmer2_orient[0]), color='black', scale=25, width=0.002)
    plt.quiver(squirmer1_x[0], squirmer1_y[0], np.cos(squirmer1_orient[0]), np.sin(squirmer1_orient[0]), color='black', scale=25, width=0.002)

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

def plot_sim_squirmer_border(R, histories, filename, dir='graphs'):
    #Plot 4 squirmers with different beta
    #histories = [[historyb01, historyb02], [historybinf1, historybinf2], [historybsup1, historybsup2]]
    plt.figure(figsize=(8, 8))
    plt.plot([-R, R], [-R, -R], 'k-', linewidth=2)
    plt.plot([-R, R], [R, R], 'k-', linewidth=2)
    plt.plot([-R, -R], [-R, R], 'k-', linewidth=2)
    plt.plot([R, R], [-R, R], 'k-', linewidth=2)

    colors = [('blue', 'cyan', 'darkcyan', 'steelblue'), ('orange', 'gold', 'goldenrod', 'orangered'), ('green', 'lime', 'limegreen', 'olive')]
    labels = ['beta=0', 'beta<0', 'beta>0']

    for history_pair, (color1, color2, color3, color4), label in zip(histories, colors, labels):
        squirmer1_x, squirmer1_y, squirmer1_orient = [], [], []
        squirmer2_x, squirmer2_y, squirmer2_orient = [], [], []
        squirmer3_x, squirmer3_y, squirmer3_orient = [], [], []
        squirmer4_x, squirmer4_y, squirmer4_orient = [], [], []

        for step in history_pair[0]:
            squirmer1_x.append(step[0])
            squirmer1_y.append(step[1])
            squirmer1_orient.append(step[4])
            squirmer2_x.append(step[2])
            squirmer2_y.append(step[3])
            squirmer2_orient.append(step[5])

        for step in history_pair[1]:
            squirmer3_x.append(step[0])
            squirmer3_y.append(step[1])
            squirmer3_orient.append(step[4])
            squirmer4_x.append(step[2])
            squirmer4_y.append(step[3])
            squirmer4_orient.append(step[5])

        plt.plot(squirmer1_x, squirmer1_y, label=f'Squirmer1 {label}', color=color1)
        plt.plot(squirmer2_x, squirmer2_y, label=f'Squirmer2 {label}', color=color2)
        plt.plot(squirmer3_x, squirmer3_y, label=f'Squirmer3 {label}', color=color3)
        plt.plot(squirmer4_x, squirmer4_y, label=f'Squirmer4 {label}', color=color4)

        last_orient1 = squirmer1_orient[0]
        last_orient2 = squirmer2_orient[0]
        last_orient3 = squirmer3_orient[0]
        last_orient4 = squirmer4_orient[0]

        for i in range(len(squirmer1_orient)):
            if squirmer1_orient[i] != last_orient1:
                plt.quiver(squirmer1_x[i], squirmer1_y[i], np.cos(squirmer1_orient[i]), np.sin(squirmer1_orient[i]), color=color1, scale=25, width=0.005)
                last_orient1 = squirmer1_orient[i]
            if squirmer2_orient[i] != last_orient2:
                plt.quiver(squirmer2_x[i], squirmer2_y[i], np.cos(squirmer2_orient[i]), np.sin(squirmer2_orient[i]), color=color2, scale=25, width=0.005)
                last_orient2 = squirmer2_orient[i]
            if squirmer3_orient[i] != last_orient3:
                plt.quiver(squirmer3_x[i], squirmer3_y[i], np.cos(squirmer3_orient[i]), np.sin(squirmer3_orient[i]), color=color3, scale=25, width=0.005)
                last_orient3 = squirmer3_orient[i]
            if squirmer4_orient[i] != last_orient4:
                plt.quiver(squirmer4_x[i], squirmer4_y[i], np.cos(squirmer4_orient[i]), np.sin(squirmer4_orient[i]), color=color4, scale=25, width=0.005)
                last_orient4 = squirmer4_orient[i]


    #Plot initial orientations
    plt.quiver(squirmer1_x[0], squirmer1_y[0], np.cos(squirmer1_orient[0]), np.sin(squirmer1_orient[0]), color='black')
    plt.quiver(squirmer2_x[0], squirmer2_y[0], np.cos(squirmer2_orient[0]), np.sin(squirmer2_orient[0]), color='black')
    plt.quiver(squirmer3_x[0], squirmer3_y[0], np.cos(squirmer3_orient[0]), np.sin(squirmer3_orient[0]), color='black')
    plt.quiver(squirmer4_x[0], squirmer4_y[0], np.cos(squirmer4_orient[0]), np.sin(squirmer4_orient[0]), color='black')

    plt.axis('equal')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Positions and Orientations of Squirmers')
    plt.legend(fontsize="7")
    plt.grid(True)

    if not os.path.exists(dir):
        os.makedirs(dir)
    save_path = os.path.join(dir, filename + '.png')
    plt.savefig(save_path)
    plt.close()

def create_video_from_history(history, R, filename='squirmers_simulation.mp4', dir='videos', fps=30, gif=False):
    if not os.path.exists(dir):
        os.makedirs(dir)
    if gif:
        filename = filename if filename.endswith('.gif') else filename + '.gif'
    else:
        filename = filename if filename.endswith('.mp4') else filename + '.mp4'
    save_path = os.path.join(dir, filename)

    fig, ax = plt.subplots()
    ax.set_xlim(-R, R)
    ax.set_ylim(-R, R)
    line1, = ax.plot([], [], 'bo', markersize=8, label='Squirmer 1')
    line2, = ax.plot([], [], 'ro', markersize=8, label='Squirmer 2')
    orientation1, = ax.plot([], [], 'b-', lw=1)
    orientation2, = ax.plot([], [], 'r-', lw=1)
    ax.legend()

    def init():
        line1.set_data([], [])
        line2.set_data([], [])
        orientation1.set_data([], [])
        orientation2.set_data([], [])
        return line1, line2, orientation1, orientation2

    def update(frame):
        x1, y1, x2, y2, theta1, theta2, *_ = history[frame]
        line1.set_data(x1, y1)
        line2.set_data(x2, y2)
        orientation1.set_data([x1, x1 + 0.1 * np.cos(theta1)], [y1, y1 + 0.1 * np.sin(theta1)])
        orientation2.set_data([x2, x2 + 0.1 * np.cos(theta2)], [y2, y2 + 0.1 * np.sin(theta2)])
        return line1, line2, orientation1, orientation2

    ani = FuncAnimation(fig, update, frames=len(history), init_func=init, blit=True)

    if gif == True:
        #For saving animation as GIF
        ani.save('squirmers_simulation.gif', writer=PillowWriter(fps=fps))
    else:
        #Save animation as video
        ani.save(save_path, writer='ffmpeg', fps=fps)