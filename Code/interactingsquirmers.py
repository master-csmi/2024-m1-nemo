import numpy as np
import matplotlib
matplotlib.use('Agg')
from squirmer import Squirmer
from plot import plot_sim_nsquirmers, create_video_from_history, plot_time

class InteractingSquirmers:

    def __init__(self, N, xs, ys, orientations, radius, beta, v0, Nx, Ny, dt, dt_out, T, Es, ds, mu, R, lnEps_cr, D, n, no, border=True):
        self.N = N
        self.xs = np.array(xs, dtype=float)
        self.ys = np.array(ys, dtype=float)
        self.orientations = np.array(orientations, dtype=float)
        self.radius = radius
        self.beta = beta
        self.v0 = v0
        self.Nx = Nx
        self.Ny = Ny
        self.dt = dt
        self.dt_out = dt_out
        self.Es = Es
        self.mu = mu
        self.R = R
        self.ds = ds
        self.T = T
        self.lnEps_cr = lnEps_cr
        self.D = D
        self.Do = 3*D/(4*radius**2)
        self.n = n
        self.no = no
        self.nos = np.zeros(N, dtype=float)
        self.Fs_x, self.Fs_y, self.Fl_x, self.Fl_y, self.val, self.gamma_w = np.zeros(N, dtype=float), np.zeros(N, dtype=float), np.zeros(N, dtype=float), np.zeros(N, dtype=float), np.zeros(N, dtype=float), np.zeros(N, dtype=float)
        self.Fs_pw = np.zeros((2,N), dtype=float)
        #border = true || false, true for reflective, false for periodic 
        self.border = border

        colors = list(matplotlib.colors.CSS4_COLORS.keys())
        colors = [color for color in colors if not self.is_light_color(matplotlib.colors.CSS4_COLORS[color])]
        np.random.shuffle(colors)
        self.colors = colors[:N]
        
        self.squirmers = np.empty(N, dtype=object)

        for i in range(N):
            self.squirmers[i] = Squirmer(self.xs[i], self.ys[i], self.orientations[i], radius, beta, v0)
        self.B1 = self.squirmers[0].B1
        self.B2 = self.squirmers[0].B2

    def is_light_color(self, hex_color):
        #Define what a color too bright is
        rgb = matplotlib.colors.hex2color(hex_color)
        luminance = 0.2126 * rgb[0] + 0.7152 * rgb[1] + 0.0722 * rgb[2]
        return luminance > 0.7

    def polar_order_parameter(self):
        """Uses the orientations of all of
        the squirmers to return the polar order parameter"""
        sum1 = 1./self.N*np.sum(np.cos(self.orientations))
        sum2 =  1./self.N*np.sum(np.sin(self.orientations))
        summ_final = np.sqrt(sum1**2 + sum2**2)
        return summ_final
    
    def is_in_square(self):
        #return True if the squirmers are in the square
        return np.all(abs(self.xs) <= (self.Nx - self.radius)) & np.all(abs(self.ys) <= (self.Ny - self.radius))

    def check_squirmers_square(self):
        #return True if the squirmers have been successfully initialized
        if (self.is_in_square() == False):
            raise ValueError("Squirmers must be inside the square")
        return True
    
    def distance_all(self):
        """Returns the distance (xj - xi) between all of the squirmers"""
        Dxs = self.xs[:, None] - self.xs
        Dys = self.ys[:, None] - self.ys
        dists = np.sqrt(Dxs**2 + Dys**2)
        return Dxs, Dys, dists

    def forcesSteric(self, Dxs, Dys, dists):
        """Uses the distances, mu, a and Es to compute and return
        Fs_x : the steric forces in the x axis
        Fs_y : the steric forces in the y axis"""
        a = self.radius

        Vc = np.minimum(a/dists, 0.5)

        tmp = -(self.Es/(np.pi*self.mu*a**2))*(2*(2*Vc)**13-(2*Vc)**7)/dists
        Fs_x = tmp*Dxs
        Fs_y = tmp*Dys
        return Fs_x, Fs_y
    
    def torquesLubrification(self, Dx, Dy, dist, theta):
        """Uses the distances, the orientations, B1, B2, a, mu and lnEps_cr
        to compute and return the lubrification torques
        val are the torques exerted on the ith squirmer by its own flow
        val2 are the torques exerted on the jth squirmer by the flow of the ith one"""
        B1 = self.B1
        B2 = self.B2
        a = self.radius

        eieijt = (np.cos(theta)*Dy - np.sin(theta)*Dx)/dist
        cosalpha = (np.cos(theta)*Dx + np.sin(theta)*Dy)/dist

        sinalpha = np.sqrt(np.maximum((1 - cosalpha * cosalpha), 0))
        somme = B1 * sinalpha + B2 * cosalpha*sinalpha

        lnEps = np.log(np.maximum(self.lnEps_cr,(dist/a - 2)))
                
        val = (8/5)*self.mu*np.pi*(a**2)*eieijt*somme*lnEps
        val2 = (1/4)*val
        
        return val, val2
        
    def forcesLubrification(self, Dx, Dy, dist, theta):
        """Uses the distances, orientations, B1, B2, mu, and a to compute and return
        F_x : the lubrification forces in the x axis
        F_y : the lubrification forces in the y axis"""
        B1 = self.B1
        B2 = self.B2
        a = self.radius

        eieijt = (np.cos(theta)*Dy - np.sin(theta)*Dx)/dist
        cosalpha = (np.cos(theta)*Dx + np.sin(theta)*Dy)/dist

        sinalpha = np.sqrt(np.maximum((1 - cosalpha * cosalpha), 0))
        somme = B1 * sinalpha + B2 * cosalpha*sinalpha
        sommeFz = B1 * sinalpha * cosalpha + (1/2)*B1 * cosalpha * eieijt**2 + B2 * sinalpha * cosalpha**2 + (1/2)*B2 * (2*cosalpha**2-1) * eieijt**2

        lnEps = np.log(np.maximum(self.lnEps_cr,(dist/a - 2)))
        
        #lambda=1
        F_x = -np.pi * self.mu * a * eieijt * somme * lnEps * Dx
        F_y = -9 * self.mu * np.pi*a*(1/4)*sommeFz* lnEps * Dy

        return F_x, F_y
    
    def compute_force_squirmer_border_x(self, xs, ys):
        """Uses the distances between the squirmers and the wall, a and mu
        to compute and return the steric forces exerted on the squirmers in the x axis"""
        RRi = np.sqrt((xs - self.Nx)**2 + (ys - self.Ny)**2)
        tmp = -((self.Es*(self.Nx - xs))/(np.pi*self.mu*self.radius**2*RRi))*(2*(self.radius/RRi)**13-(self.radius/RRi)**7)
        return tmp*xs
    
    def compute_force_squirmer_border_y(self, xs, ys):
        """Uses the distances between the squirmers and the wall, a and mu
        to compute and return the steric forces exerted on the squirmers in the y axis"""
        RRi = np.sqrt((xs - self.Nx)**2 + (ys - self.Ny)**2)
        tmp = -((self.Es*(self.Ny - ys))/(np.pi*self.mu*self.radius**2*RRi))*(2*(self.radius/RRi)**13-(self.radius/RRi)**7)
        return tmp*ys
    
    def force_torque_lubrification_border_x(self, xs, theta, border):
        """Uses the distances between the squirmers and the wall, a, B1, B2 and mu
        to compute and return the lubrification forces and torques
        exerted on the squirmers in the x axis"""
        B1 = self.B1
        B2 = self.B2
        a = self.radius

        if border == 1:
            #border = 1 for right border
            Dx = self.Nx - xs
        else:
            #else for left border
            Dx = -self.Nx - xs
        Dy = 0
        dist = np.sqrt(Dx**2 + Dy**2)

        eieijt = (np.cos(theta)*Dy - np.sin(theta)*Dx)/dist
        cosalpha = (np.cos(theta)*Dx + np.sin(theta)*Dy)/dist

        sinalpha = np.sqrt(np.maximum((1 - cosalpha * cosalpha), 0))
        somme = B1 * sinalpha + B2 * cosalpha*sinalpha
        lnEps = np.log(np.maximum(self.lnEps_cr,(dist/a - 2)))
        
        #lambda=1
        F_x = -(4/5) * np.pi * self.mu * a * eieijt * somme * lnEps * Dx
        gamma_w = (16/5)*self.mu*np.pi*(a**2)*eieijt*somme*lnEps
        return F_x, gamma_w
    
    def force_torque_lubrification_border_y(self, ys, theta, border):
        """Uses the distances between the squirmers and the wall, a, B1, B2 and mu
        to compute and return the lubrification forces and torques
        exerted on the squirmers in the y axis"""
        B1 = self.B1
        B2 = self.B2
        a = self.radius
        if border == 1:
            #border = 1 for upper wall
            Dy = self.Ny - ys
        else:
            #else for lower wall
            Dy = -self.Ny - ys
        Dx = 0
        dist = np.sqrt(Dx**2 + Dy**2)

        eieijt = (np.cos(theta)*Dy - np.sin(theta)*Dx)/dist
        cosalpha = (np.cos(theta)*Dx + np.sin(theta)*Dy)/dist

        sinalpha = np.sqrt(np.maximum((1 - cosalpha * cosalpha), 0))
        somme = B1 * sinalpha + B2 * cosalpha*sinalpha
        sommeFz = B1 * sinalpha * cosalpha + (1/2)*B1 * cosalpha * eieijt**2 + B2 * sinalpha * cosalpha**2 + (1/2)*B2 * (2*cosalpha**2-1) * eieijt**2
        lnEps = np.log(np.maximum(self.lnEps_cr,(dist/a - 2)))
        
        #lambda=1
        F_y = -9 * self.mu * np.pi*a*sommeFz* lnEps * Dy
        gamma_w = (16/5)*self.mu*np.pi*(a**2)*eieijt*somme*lnEps
        return F_y, gamma_w

    #Reflective boundary condition
    def ref_border_x(self, xs, orientation, boundary):
        """takes the distance between the squirmers and the x-wall
        and returns x-wall - distance"""
        orientation = np.pi - orientation
        if boundary == 1:
            #1 for the right border
            diff = xs + self.radius - self.Nx
            xs = self.Nx - diff - self.radius
        else:
            diff = xs - self.radius + self.Nx
            xs = -self.Nx - diff + self.radius

        return xs, orientation
    
    def ref_border_y(self, ys, orientation, boundary):
        """takes the distance between the squirmers and the y-wall
        and returns y-wall - distance"""
        orientation = -orientation
        if boundary == 1:
            #1 for the up boundary
            diff = ys + self.radius - self.Ny
            ys = self.Ny - diff - self.radius
        else:
            diff = ys - self.radius + self.Ny
            ys = -self.Ny - diff + self.radius

        return ys, orientation

    def perio_border_x(self, xs, boundary):
        """returns the xs that simulates a periodic border"""
        if boundary == 1:
            #1 for right boundary
            xs -= 2*self.Nx
        else:
            xs += 2*self.Nx
        return xs

    def loop_time(self):
        """simulates the behaviors of the squirmers
        for T final time with a time-step dt and puts the informations
        in a list with a time-step dt_out, it returns the list"""
        self.vector_dists_min = []
        tout = self.dt_out
        a = self.radius
        history = []
        data = [self.xs.tolist(), self.ys.tolist(), self.orientations.tolist(),
                self.Fs_x.tolist(), self.Fs_y.tolist(), self.Fl_x.tolist(), self.Fl_y.tolist(),
                self.val.tolist(), self.gamma_w.tolist(), self.Fs_pw.tolist(), 0]
        history.append(data)
        self.list_polar = []
        self.list_cluster_param = []

        for t in np.arange(0, self.T, self.dt):
            self.Fs_x.fill(0)
            self.Fs_y.fill(0)
            self.Fl_x.fill(0)
            self.Fl_y.fill(0)
            self.val.fill(0)
            self.gamma_w.fill(0)
            self.Fs_pw.fill(0)
            self.nos.fill(0)

            Dxs, Dys, dists = self.distance_all()

            #Minimum distance between all of the squirmers
            dist = np.array(dists)
            dist_nz = dist[dist!=0]
            if dist_nz.size > 0:
                min_dist = np.min(dist_nz - 2*a)
                self.vector_dists_min.append(min_dist)
                # if min_dist < 0:
                    # print(f"min_dist = {min_dist}")
                    # print(f"t_min_dist <0 = {t}")
                    # indices_min = np.argwhere((dist - 2 * a) == min_dist)

                    # print(f"Squirmer1: x={self.xs[indices_min[0][0]]}, y={self.ys[indices_min[0][1]]}")
                    # print(f"Squirmer2: x={self.xs[indices_min[1][0]]}, y={self.ys[indices_min[1][1]]}")
            
            #Clustering order parameter
            dist_neigh = (dists<self.R)&(dists!=0)
            n_neigh = np.sum(dist_neigh, axis=1)
            self.list_cluster_param.append((1/(self.N*self.N/2))*sum(n_neigh))

            #Steric Forces
            dist_steric = (dists<self.ds)&(dists!=0)
            dist_lubrification = (dists<=3*a)&(dists!=0)
            j_dist_steric = np.where(dist_steric)
            j_dist_lubr = np.where(dist_lubrification)
            # print(f"Dxs = {Dxs}\n Dys = {Dys}\n dists = {dists}")
            Fs_x, Fs_y = self.forcesSteric(Dxs[dist_steric], Dys[dist_steric], dists[dist_steric])
            # print(f"Fs_x = {Fs_x} and shape = {Fs_x.shape}")
            # print(f"self_Fs_x = {self.Fs_x[j_dist_steric[0]]}")
            self.Fs_x[j_dist_steric[0]] += (Fs_x)
            # print(f"self.Fs_x après = {self.Fs_x[j_dist_steric[0]]} and shape = {self.Fs_x[j_dist_steric[0]].shape}")
            # print("\n")
            self.Fs_y[j_dist_steric[0]] += (Fs_y)

            #Lubrification Forces and Torques
            Fl_x, Fl_y = self.forcesLubrification(Dxs[j_dist_lubr], Dys[j_dist_lubr], dists[j_dist_lubr], self.orientations[j_dist_lubr[0]])
            val1, val2 = self.torquesLubrification(Dxs[j_dist_lubr], Dys[j_dist_lubr], dists[j_dist_lubr], self.orientations[j_dist_lubr[0]])
            self.Fl_x[j_dist_lubr[0]] += Fl_x
            self.Fl_y[j_dist_lubr[0]] += Fl_y
            self.Fl_x[j_dist_lubr[1]] -= Fl_x
            self.Fl_y[j_dist_lubr[1]] -= Fl_y
            self.val[j_dist_lubr[0]] += val1
            self.val[j_dist_lubr[1]] += val2
            # print(f"Fl_x = {self.Fl_x}")
            # print(f"Fl_y = {self.Fl_y}\n")

            #Noises
            self.ns = np.random.uniform(-self.n/2, self.n/2, size=self.N)
            self.nos = np.random.uniform(-self.no/2, self.no/2, size=self.N)

            #Forces and torques between a squirmer and a border
            dist_forces_x = (self.Nx-abs(self.xs)) < 2**(1/6)*a
            dist_forces_y = (self.Ny-abs(self.ys)) < 2**(1/6)*a

            dist_lub_x = (self.Nx-abs(self.xs)) < 2*a
            dist_lub_y = (self.Ny-abs(self.ys)) < 2*a

            #If we simulate in a box
            if self.border:
                i_dist_force_x = np.where(dist_forces_x)[0]
                i_dist_lub_x = np.where(dist_lub_x)[0]
                
                self.Fs_pw[0][i_dist_force_x] += self.compute_force_squirmer_border_x(self.xs[i_dist_force_x], self.ys[i_dist_force_x])

                i_dist_lub_x_r = i_dist_lub_x[self.xs[i_dist_lub_x] > 0]
                i_dist_lub_x_l = i_dist_lub_x[self.xs[i_dist_lub_x] < 0]
                
                Fl_xr, gamma_wr = self.force_torque_lubrification_border_x(self.xs[i_dist_lub_x_r], self.orientations[i_dist_lub_x_r], 1)
                Fl_xl, gamma_wl = self.force_torque_lubrification_border_x(self.xs[i_dist_lub_x_l], self.orientations[i_dist_lub_x_l], 2)
                
                self.Fl_x[i_dist_lub_x_r] += Fl_xr
                self.Fl_x[i_dist_lub_x_l] += Fl_xl
                self.gamma_w[i_dist_lub_x_r] += gamma_wr
                self.gamma_w[i_dist_lub_x_l] += gamma_wl

            i_dist_force_y = np.where(dist_forces_y)[0]
            i_dist_lub_y = np.where(dist_lub_y)[0]

            self.Fs_pw[1][i_dist_force_y] += self.compute_force_squirmer_border_y(self.xs[i_dist_force_y], self.ys[i_dist_force_y])
            
            i_dist_lub_y_u = i_dist_lub_y[self.ys[i_dist_lub_y] > 0]
            i_dist_lub_y_d = i_dist_lub_y[self.ys[i_dist_lub_y] < 0]
            
            Fl_yu, gamma_wu = self.force_torque_lubrification_border_y(self.ys[i_dist_lub_y_u], self.orientations[i_dist_lub_y_u], 1)
            Fl_yd, gamma_wd = self.force_torque_lubrification_border_y(self.ys[i_dist_lub_y_d], self.orientations[i_dist_lub_y_d], 2)
            
            self.Fl_y[i_dist_lub_y_u] += Fl_yu
            self.Fl_y[i_dist_lub_y_d] += Fl_yd
            self.gamma_w[i_dist_lub_y_u] += gamma_wu
            self.gamma_w[i_dist_lub_y_d] += gamma_wd

            # print(self.val)

            #Evolution of position
            self.orientations += self.dt*(self.val + self.gamma_w) + np.sqrt(2*self.dt*self.Do)*self.nos
            self.xs += self.dt*(self.v0*np.cos(self.orientations) - self.Fs_x - self.Fs_pw[0] + self.Fl_x) + np.sqrt(2*self.dt*self.D)*self.ns
            self.ys += self.dt*(self.v0*np.sin(self.orientations) - self.Fs_y - self.Fs_pw[1] + self.Fl_y) + np.sqrt(2*self.dt*self.D)*self.ns

            #Borders
            mask_x1 = (self.xs + a) > self.Nx
            mask_x2 = (self.xs - a) < -self.Nx
            mask_y1 = (self.ys + a) > self.Ny
            mask_y2 = (self.ys - a) < -self.Ny
            #x_borders
            if np.any(mask_x1):
                #Box simulation
                if self.border:
                    self.xs[mask_x1], self.orientations[mask_x1] = self.ref_border_x(self.xs[mask_x1], self.orientations[mask_x1], 1)
                #Chanel simulation
                else:
                    self.xs[mask_x1] = self.perio_border_x(self.xs[mask_x1], 1)
            if np.any(mask_x2):
                #Box simulation
                if self.border:
                    self.xs[mask_x2], self.orientations[mask_x2] = self.ref_border_x(self.xs[mask_x2], self.orientations[mask_x2], 2)
                #Chanel simulation
                else:
                    self.xs[mask_x2] = self.perio_border_x(self.xs[mask_x2], 2)
            
            #y_borders
            if np.any(mask_y1):
                self.ys[mask_y1], self.orientations[mask_y1] = self.ref_border_y(self.ys[mask_y1], self.orientations[mask_y1], 1)
            if np.any(mask_y2):
                self.ys[mask_y2], self.orientations[mask_y2] = self.ref_border_y(self.ys[mask_y2], self.orientations[mask_y2], 2)

            for i, s in enumerate(self.squirmers):
                s.x = self.xs[i]
                s.y = self.ys[i]
                s.orientation = self.orientations[i]
            
            #Polar order parameter
            self.list_polar.append(self.polar_order_parameter())
            
            if t >= tout:
                print(f"{(tout/self.T)*100} %")
                data = [self.xs.tolist(), self.ys.tolist(), self.orientations.tolist(),
                        self.Fs_x.tolist(), self.Fs_y.tolist(), self.Fl_x.tolist(), self.Fl_y.tolist(),
                        self.val.tolist(), self.gamma_w.tolist(), self.Fs_pw.tolist(), tout]
                history.append(data)
                tout += self.dt_out

        self.history = history
        return history

def run(choice, N, a, beta, v0, Nx, Ny, dt, dt_out, T, Es, ds, mu, R, lnEps_cr, D, n, no, border, filename, border_plot=False):
    """
    Runs the simulations depending on the 'choice' parameter
    choice: 'plot' to plot the graph of the behaviors of the squirmers in the 'graphs' directory
            'video' to create a video of the behaviors of the squirmers in the 'videos' directory
            'Eo_sim', not used anymore, to test every Eo parameter and plot or make a video
            'border', to simulate the behaviors of a squirmer near a wall, only plot possible
    N: number of squirmers
    a: radius of the squirmers
    beta: beta of the squirmers
    v0: velocity of the squirmers
    Nx, Ny: length and height of the rectangle
    dt, dt_out: time-step of simulation, time-step of information
    T: final time of simulation
    Es: amplitude of steric forces
    ds: distance of steric interactions
    mu: dynamic viscosity
    R: Distance where other squirmers seen as neighbour
    lnEps_cr: cut-off for -log
    D: Translational diffusivity
    n: Translational noise
    no: angular noise
    border: False to simulate a chanel and True to simulate a box
    filename: used only with 'plot' and 'video' choices, name of the file
              in which the videos or plots will be saved
    border_plot: used only with 'plot', False to not plot the borders and True to plot them
    """
    if choice == 'video':
        #coordinates and orientations
        orients = np.zeros(N, dtype=float)
        orients = np.random.uniform(0, 2*np.pi, size=N)
        xs = np.empty(N)
        ys = np.empty(N)
        dir = 'videos'

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
        interact = InteractingSquirmers(N, xs, ys, orients, a, beta, v0, Nx, Ny, dt, dt_out, T, Es, ds, mu, R, lnEps_cr, D, n, no, border)
        history = interact.loop_time()
        plot_time(interact, interact.vector_dists_min, "min_dist_" + filename, 'minimal distance', dir)
        plot_time(interact, interact.list_polar, "polar_" + filename, 'polar parameter', dir)
        plot_time(interact, interact.list_cluster_param, "cluster_" + filename, 'clustering order parameter', dir)
        print(f"Mean polar order = {np.mean(interact.list_polar)}")

        create_video_from_history(history, Nx, Ny, N, a, filename=filename, dir=dir)
    elif choice == 'plot':
        #coordinates and orientations
        orients = np.zeros(N, dtype=float)
        orients = np.random.uniform(0, 2*np.pi, size=N)
        xs = np.empty(N)
        ys = np.empty(N)
        dir = 'graphs'

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
        interact = InteractingSquirmers(N, xs, ys, orients, a, beta, v0, Nx, Ny, dt, dt_out, T, Es, ds, mu, R, lnEps_cr, D, n, no, border)
        history = interact.loop_time()
        plot_time(interact, interact.vector_dists_min, "min_dist_" + filename, 'minimal distance', dir=dir)
        plot_time(interact, interact.list_polar, "polar_" + filename, 'polar parameter', dir=dir)
        plot_time(interact, interact.list_cluster_param, "cluster_" + filename, 'clustering order parameter', dir=dir)

        plot_sim_nsquirmers(history, Nx, Ny, N, a, border_plot, False, filename=filename, dir=dir)
    elif choice == 'border':
        a = 0.05
        xs = [-0.4]
        ys = [-0.7]
        orient = [[-np.pi/6], [-np.pi/4], [-np.pi/3], [-np.pi/2]] 
        N = 1
        sim_border = True
        T = 0.9
        v0 = 1
        Nx = 0.5
        Ny = 1
        for i, pi in enumerate(orient):
            filename = 'sim_num_' + str(i)
            print(filename)
            interact = InteractingSquirmers(N, xs, ys, pi, a, beta, v0, Nx, Ny, dt, dt_out, T, Es, ds, mu, R, lnEps_cr, D, n, no, border)
            history = interact.loop_time()

            plot_sim_nsquirmers(history, Nx, Ny, N, a, border_plot, sim_border, filename=filename, dir='graphs/border')
    elif choice == 'Eo_sim':
        #amplitude of orientational interactions
        #Eo[0] = Eoinitial, E[1] = Eobrumley, E[2] = Eolauga
        Eo = [(((3./10.)*v0/a), "Eo_init"), ((16/10)*mu*np.pi*a**2, "Eo_brumley"), ((-3./2.)*(v0/a), "Eo_lauga"),
            (((-3./10.)*v0/a), "mEo_init"), (-5, "m5"), (0.005, "m0_005"), (-2, "m2"), (-0.5, "m0_5"),
            (0.5, "0_5")]
        betas = [(0, "beta0"), (-7.5, "betainf"), (7.5, "betasup")]
        while True:
            output_type = input("Which type of simulation? (plot/video): ").strip().lower()
            if output_type in ['plot', 'video']:
                break
            else:
                print("Invalid input. Please enter 'plot' or 'video'.")
        xseo = [-a, 2*a/1.5]
        yseo = [0, 0]
        orient1 = np.pi/2
        orient2 = [(np.pi/2, "pi_2_"), (-np.pi/2, "mpi_2_"), (3*np.pi/4, "3pi_4_"), (-3*np.pi/4, "m3pi_4_"), (np.pi, "pi_"), (2*np.pi, "2pi_"), (np.pi/4, "pi_4_"), (-np.pi/4, "mpi_4_")]
        for (betaeo, labelbeta) in betas:
            for (Eos, labeleo) in Eo:
                for (pi, labelpi) in orient2:
                    filenameeo = labelpi + labeleo
                    orientseo = [orient1, pi]
                    interact = InteractingSquirmers(N, xseo, yseo, orientseo, a, betaeo, v0, Nx, Ny, dt, dt_out, T, Es, ds, mu, Eos, lnEps_cr, D, n, no, border)
                    history = interact.loop_time()
                    if output_type == 'plot':
                        direo = 'graphs/Eo_analysis/' + labelbeta + '/' + labelpi
                        plot_sim_nsquirmers(history, 1, 1, 2, a, border_plot, False, filename=filenameeo, dir=direo)
                    else:
                        direo = 'videos/Eo_analysis/' + labelbeta + '/' + labelpi
                        create_video_from_history(history, 1, 1, 2, a, filename=filenameeo, dir=direo)