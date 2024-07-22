from Code.interactingsquirmers import InteractingSquirmers
import numpy as np

def test_interact():
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
    a = 0.05
    #betas
    beta = 0
    #time-step
    dt = 1e-4
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
    ds = 2*a*2**(1./6)
    #Translational diffusivity
    D = 0
    #Translational noise
    n = 1e-2
    #Angular noise
    no = 1e-2
    #Distance of particle seen as "Neighbour"
    R = 0.07
    #border defines the simulation, in a chanel(False) or a box(True)
    border = True

    N = 375

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

    interact = InteractingSquirmers(N, xs, ys, orients, a, beta, v0, Nx, Ny, dt, dt_out, T, Es, ds, mu, R, lnEps_cr, D, n, no, border)
    interact.loop_time()
    assert np.all(interact.vector_dists_min > 0)