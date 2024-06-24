import matplotlib
import numpy as np
import os
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def is_light_color(hex_color):
    #Define what a color too bright is
    rgb = matplotlib.colors.hex2color(hex_color)
    luminance = 0.2126 * rgb[0] + 0.7152 * rgb[1] + 0.0722 * rgb[2]
    return luminance > 0.7

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

def plot_sim_nsquirmers(histories, R, N, a, border, filename, dir='graphs'):
    #a : radius of the squirmer
    fig = plt.figure(figsize=(8, 8))

    if border == True:
        plt.plot([-R, R], [-R, -R], 'k-', linewidth=2)
        plt.plot([-R, R], [R, R], 'k-', linewidth=2)
        plt.plot([-R, -R], [-R, R], 'k-', linewidth=2)
        plt.plot([R, R], [-R, R], 'k-', linewidth=2)

    colors = ['blue', 'cyan', 'orange', 'gold', 'green', 'lime', 'red', 'pink', 'purple', 'violet']

    xs = [history[0] for history in histories]
    ys = [history[1] for history in histories]
    orientations = [history[2] for history in histories]
    squirmer_xs = []
    squirmer_ys = []
    squirmer_orients = []
    for i in range(N):
        x = [step[i] for step in xs]
        y = [step[i] for step in ys]
        orient = [step[i] for step in orientations]
        squirmer_xs.append(x)
        squirmer_ys.append(y)
        squirmer_orients.append(orient)

    nb_pixl_fig = fig.get_size_inches()[0]*fig.dpi
    radius_scatter = nb_pixl_fig/(2*R/a)
    s = radius_scatter**2
    scale_arrow = 15
    w = 0.007

    for i in range(N):
        plt.plot(squirmer_xs[i], squirmer_ys[i], color=colors[i])
        last_orient = squirmer_orients[i][0]
        plot_circle = 0
        for j in range(len(squirmer_orients[i])):
            new_orient = squirmer_orients[i][j]
            if new_orient != last_orient:
                plt.quiver(squirmer_xs[i][j], squirmer_ys[i][j], np.cos(new_orient), np.sin(new_orient), color=colors[i], scale=scale_arrow, width=w)
                last_orient = new_orient
                plot_circle += 1
                if plot_circle == 4:
                    plt.scatter(squirmer_xs[i][j], squirmer_ys[i][j], color=colors[i], s=s)
                    plot_circle = 0

    #Plot initial orientations
    xs = histories[0][0]
    ys = histories[0][1]
    orientations = histories[0][2]
    for i in range(N):
        plt.scatter(xs[i], ys[i], color=colors[i], s=s)
        plt.quiver(xs[i], ys[i], np.cos(orientations[i]), np.sin(orientations[i]), color='black', scale=scale_arrow, width=w)

    plt.scatter([-R, R], [-R, R], color='white', alpha=0)

    plt.axis('equal')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Positions and Orientations of Squirmers')
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
                plt.scatter(squirmer1_x[i], squirmer1_y[i], color=color1)
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

def create_video_from_history(history, R, N, filename='squirmers_simulation.mp4', dir='videos', fps=30):
    if not os.path.exists(dir):
        os.makedirs(dir)
    filename = filename if filename.endswith('.mp4') else filename + '.mp4'
    save_path = os.path.join(dir, filename)

    fig, ax = plt.subplots()
    ax.set_xlim(-R, R)
    ax.set_ylim(-R, R)
    
    lines = []
    orientations = []
    colors = ['b', 'r', 'g', 'c', 'm', 'y', 'k']
    for i in range(N):
        color = colors[i % len(colors)]
        line, = ax.plot([], [], color + 'o', markersize=8, label=f'Squirmer {i + 1}')
        orientation, = ax.plot([], [], color + '-', lw=1)
        lines.append(line)
        orientations.append(orientation)
    
    ax.legend()

    def init():
        for line in lines:
            line.set_data([], [])
        for orientation in orientations:
            orientation.set_data([], [])
        return lines + orientations

    def update(frame):
        xs, ys, thetas = history[frame][:3]
        for i in range(N):
            lines[i].set_data(xs[i], ys[i])
            orientations[i].set_data(
                [xs[i], xs[i] + 0.1 * np.cos(thetas[i])], 
                [ys[i], ys[i] + 0.1 * np.sin(thetas[i])]
            )
        return lines + orientations

    ani = FuncAnimation(fig, update, frames=len(history), init_func=init, blit=True)
    # Save animation as video
    ani.save(save_path, writer='ffmpeg', fps=fps)
    plt.close()