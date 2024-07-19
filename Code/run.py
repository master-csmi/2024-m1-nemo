import argparse
import numpy as np
import time
from interactingsquirmers import run
from simulation import sim_vicsek

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
    Es = 0.1
    #simulation time
    T = 2
    #periodicity of outputs
    dt_out = 0.01
    #viscovity parameter
    mu = 0.1
    #distance of steric interactions
    ds = 2**(7./6)*a
    #Translational diffusivity
    D = 0
    #angular noise
    no = 1e-2
    #Distance of particle seen as "Neighbour"
    R = 0.07

    #border to simulate chanel or box
    border = False
    #border_plot to plot the borders when plotting or not
    border_plot = False

    if simulation == 'video':
        run('video', N, a, beta, v0, Nx, Ny, dt, dt_out, T, Es, ds, mu, R, lnEps_cr, D, no, border, filename, border_plot)
    elif simulation == 'plot':
        run('plot', N, a, beta, v0, Nx, Ny, dt, dt_out, T, Es, ds, mu, R, lnEps_cr, D, no, border, filename, border_plot)
    elif simulation == 'border':
        run('border', N, a, beta, v0, Nx, Ny, dt, dt_out, T, Es, ds, mu, R, lnEps_cr, D, no, border, filename, border_plot)
    elif simulation == 'Eo_sim':
        run('Eo_sim', N, a, beta, v0, Nx, Ny, dt, dt_out, T, Es, ds, mu, R, lnEps_cr, D, no, border, filename, border_plot)
    elif simulation == 'vicsek':
        sim_vicsek()
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
