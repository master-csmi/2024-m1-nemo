import argparse
import numpy as np
import time
from simulation import sim_interacting_squirmers, sim_vid_interact_sq, sim_Eo_param, sim_vicsek

def main(simulation, N, filename):
    start_time = time.time()
    # Define parameters
    #velocity #F*dt < v0
    v0 = 0.3
    #length x and y axis
    lbda = 2
    Ly = 3
    Lx = lbda*Ly
    #half of the length of axis
    Nx = Lx/2
    Ny = Ly/2
    #squirmers' radius
    a = 0.02
    #betas
    beta = 0
    #time-step
    dt = 1e-3
    #cut-off for -log
    lnEps_cr = a*0.001
    #amplitude of steric interactions
    Es = 0.5
    #simulation time
    T = 100
    #periodicity of outputs
    dt_out = 0.05
    #viscovity parameter
    mu = 1
    #amplitude of orientational interactions
    #Eo[0] = Eoinitial, E[1] = Eobrumley, E[2] = Eolauga
    Eo = [(((3./10.)*v0/a), "Eo_init"), ((16/10)*mu*np.pi*a**2, "Eo_brumley"), ((-3./2.)*(v0/a), "Eo_lauga"),
          (((-3./10.)*v0/a), "mEo_init"), (-5, "m5"), (0.005, "m0_005"), (-2, "m2"), (-0.5, "m0_5"),
          (0.5, "0_5")]
    #distance of steric interactions
    ds = 2**(7./6)*a
    #angular diffusivity
    Do = 0.1
    #angular noise
    no = 1e-2
    #Distance of particle seen as "Neighbour"
    R = 0.07

    #coordinates and orientations
    orients = np.zeros(N, dtype=float)
    orients = np.random.uniform(0, 2*np.pi, size=N)
    xs = np.empty(N)
    ys = np.empty(N)

    for k in range(N):
        while True:
            x = np.random.uniform(-(Nx-2*a), (Nx-2*a))
            y = np.random.uniform(-(Ny-2*a), (Ny-2*a))
            if k == 0 or np.all(np.sqrt((xs[:k] - x)**2 + (ys[:k] - y)**2) > 2*a):
                xs[k] = x
                ys[k] = y
                break
    print(xs)
    print("initialisation done")
    #border to simulate chanel or box
    border = False
    #border_plot to do the simulation of the border or not
    border_plot = True

    if simulation == 'video':
        sim_vid_interact_sq(N, xs, ys, orients, a, beta, v0, Nx, Ny, dt, dt_out, T, Es, ds, mu, R, lnEps_cr, Do, no, border, filename, dir='videos')
    elif simulation == 'plot':
        sim_interacting_squirmers(N, xs, ys, orients, a, beta, v0, Nx, Ny, dt, dt_out, T, Es, ds, mu, R, lnEps_cr, Do, no, border, False, filename, border_plot, dir='graphs')
    elif simulation == 'border':
        a = 0.05
        xs = [-0.4]
        ys = [-0.7]
        orient = [[-np.pi/6], [-np.pi/4], [-np.pi/3], [-np.pi/2]] 
        N = 1
        sim_border = True
        T = 0.7
        v0 = 1
        for i, pi in enumerate(orient):
            filename = 'sim_num_' + str(i)
            sim_interacting_squirmers(N, xs, ys, pi, a, beta, v0, 0.5, 1, dt, dt_out, T, Es, ds, mu, R, lnEps_cr, Do, no, border, sim_border, filename, border_plot, dir='graphs/border')
    elif simulation == 'Eo_sim':
        sim_Eo_param(Eo, a, v0, dt, dt_out, T, Es, ds, mu, lnEps_cr, Do, no, border, border_plot)
    elif simulation == 'vicsek':
        N = 20
        R = 0.25
        L = 10.0
        v0 = 1.0
        beta = 0.5
        radius = 0.1
        T = 1
        dt = 0.1
        noise = 1e-4
        nb_step = 5
        sim_vicsek(N, R, L, v0, beta, radius, T, dt, noise, nb_step)
    end_time = time.time()
    print(f"Simulation time : {end_time - start_time}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run squirmer simulations.")

    parser.add_argument('simulation', choices=['video', 'plot', 'Eo_sim','border', 'vicsek'],
                        help="Choose which simulation to run")
    parser.add_argument('N', type=int, help="Number of squirmer")
    parser.add_argument('filename', type=str, help="Filename for saving the results")
    args = parser.parse_args()

    main(args.simulation, args.N, args.filename)
