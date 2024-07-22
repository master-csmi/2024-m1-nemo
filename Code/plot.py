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

def plot_sim_nsquirmers(histories, Nx, Ny, N, a, border_plot, sim_border, filename, dir='graphs'):
    #The simulation with a border is only with one squirmer
    if sim_border:
        assert N == 1

    #a : radius of the squirmer
    fig = plt.figure(figsize=(8, 8))

    if sim_border == True:
        plt.plot([-Nx, Nx], [-Ny, -Ny], 'k-', linewidth=2)
    if (border_plot == True) and (sim_border == False):
        plt.plot([-Nx, Nx], [-Ny, -Ny], 'k-', linewidth=2)
        plt.plot([-Nx, Nx], [Ny, Ny], 'k-', linewidth=2)
        plt.plot([-Nx, -Nx], [-Ny, Ny], 'k-', linewidth=2)
        plt.plot([Nx, Nx], [-Ny, Ny], 'k-', linewidth=2)

    colors = ['blue', 'cyan', 'orange', 'gold', 'green', 'lime', 'red', 'pink', 'purple', 'violet']

    xs = [history[0] for history in histories]
    ys = [history[1] for history in histories]
    orientations = [history[2] for history in histories]

    if sim_border:
        initial_position = ys[0][0]
        time = [history[-1] for history in histories]

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

    nb_pixl_fig = fig.get_size_inches()[1]*fig.dpi
    radius_scatter = nb_pixl_fig/(2*Ny/a)
    if sim_border == True or border_plot == False:
        s = radius_scatter**2
    else:
        s = radius_scatter
    scale_arrow = 15
    w = 0.005

    if N == 2:
        hypo_xs = [[], []]
        hypo_ys = [[], []]
        dt = histories[1][-1] - histories[0][-1]  # Assuming constant time step
        for i in range(N):
            x, y, orient = xs[0][i], ys[0][i], orientations[0][i]
            for t in range(len(histories)):
                hypo_xs[i].append(x)
                hypo_ys[i].append(y)
                x += np.cos(orient) * dt
                y += np.sin(orient) * dt

        # Plot hypothetical trajectories in grey
        for i in range(N):
            plt.plot(hypo_xs[i], hypo_ys[i], color='grey', linestyle='--')

    for i in range(N):
        plt.plot(squirmer_xs[i], squirmer_ys[i], color=colors[i % len(colors)])
        last_orient = squirmer_orients[i][0]
        plot_circle = 0
        reach_init_y = False
        for j in range(len(squirmer_orients[i])):
            new_orient = squirmer_orients[i][j]
            if new_orient != last_orient:
                plt.quiver(squirmer_xs[i][j], squirmer_ys[i][j], np.cos(new_orient), np.sin(new_orient), color=colors[i % len(colors)], scale=scale_arrow, width=w)
                last_orient = new_orient
                plot_circle += 1
                if plot_circle == 4:
                    plt.scatter(squirmer_xs[i][j], squirmer_ys[i][j], color=colors[i % len(colors)], s=s)
                    plot_circle = 0
            if j>0 and sim_border and not reach_init_y and squirmer_ys[i][j] >= initial_position:
                reach_init_y = True
                plt.scatter(squirmer_xs[i][j], squirmer_ys[i][j], color='red', s=s)
                plt.text(squirmer_xs[i][j] + 0.1, squirmer_ys[i][j], f'Time: {time[j]:.2f}', fontsize=12, color='red')

    #Plot initial orientations
    xs = histories[0][0]
    ys = histories[0][1]
    orientations = histories[0][2]
    for i in range(N):
        plt.scatter(xs[i], ys[i], color=colors[i % len(colors)], s=s)
        plt.quiver(xs[i], ys[i], np.cos(orientations[i]), np.sin(orientations[i]), color='black', scale=scale_arrow, width=w)

    if sim_border != True:
        plt.scatter([-Nx, Nx], [-Ny, Ny], color='white', alpha=0)

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

def plot_time(interact, list_plot, filename, label, dir='graphs'):
    plt.figure()
    t = np.arange(0, interact.T, interact.dt)
    plt.plot(t, list_plot)
    plt.xlabel('Time step')
    plt.ylabel(label)
    plt.grid(True)

    if not os.path.exists(dir):
        os.makedirs(dir)
    save_path = os.path.join(dir, filename + '.png')
    plt.savefig(save_path)
    plt.close()

def create_video_from_history(history, Nx, Ny, N, a, filename='squirmers_simulation.mp4', dir='videos', fps=30):
    if not os.path.exists(dir):
        os.makedirs(dir)
    filename = filename if filename.endswith('.mp4') else filename + '.mp4'
    save_path = os.path.join(dir, filename)

    fig, ax = plt.subplots()
    ax.set_xlim(-Nx, Nx)
    ax.set_ylim(-Ny, Ny)

    nb_pixl_fig = fig.get_size_inches()[0]*fig.dpi
    radius_scatter = nb_pixl_fig/(2*Ny/a)
    s = radius_scatter/2

    lines = []
    orientations = []
    colors = ['b', 'r', 'g', 'c', 'm', 'y', 'k']
    for i in range(N):
        color = colors[i % len(colors)]
        line, = ax.plot([], [], color + 'o', markersize=s)
        orientation, = ax.plot([], [], color + '-', lw=1)
        lines.append(line)
        orientations.append(orientation)

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