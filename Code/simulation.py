import numpy as np
import time
import matplotlib.pyplot as plt
from vicsek import Vicsek_continous

def sim_vicsek(N):
    R = 0.25
    L = 10.0
    v0 = 1.0
    beta = 0
    radius = 0.1
    T = 1
    dt = 0.1
    no = 1e-4
    nb_step = 5

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