import numpy as np
import matplotlib
matplotlib.use('Agg')
from squirmer import Squirmer

class InteractingSquirmers:

    def __init__(self, N, xs, ys, orientations, radius, beta, v0, Nx, Ny, dt, dt_out, T, Es, ds, mu, R, lnEps_cr, Do, no, border=True):
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
        self.Do = Do
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
        #Returns the polar order parameter
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
    
    def distance_sq(self, squirmer1, squirmer2):
        #return distance between the two squirmers
        Dx = squirmer2.x - squirmer1.x
        Dy = squirmer2.y - squirmer1.y
        return Dx, Dy, np.sqrt(Dx**2 + Dy**2)
    
    def distance_all(self):
        Dxs = self.xs[:, None] - self.xs
        Dys = self.ys[:, None] - self.ys
        dists = np.sqrt(Dxs**2 + Dys**2)
        return Dxs, Dys, dists
    
    def distance_center(self, squirmer):
        #Compute the distance between the squirmers and the center (0,0)
        dist = np.sqrt(squirmer.x**2 + squirmer.y**2)
        return dist

    def forcesSteric(self, Dxs, Dys, dists):
        #Compute the steric forces between two particles
        a = self.radius

        Vc = np.minimum(a/dists, 0.5)

        tmp = -3*(self.Es/(6*np.pi*self.mu*a))*(2*(2*Vc)**13-(2*Vc)**7)/np.sqrt(dists)
        Fs_x = tmp*Dxs
        Fs_y = tmp*Dys
        return Fs_x, Fs_y
    
    def torquesLubrification(self, Dx, Dy, dist, theta):
        #Computes the lubrification torques produced by two interacting squirmers
        B1 = self.B1
        B2 = self.B2
        a = self.radius

        eieijt = (np.cos(theta)*Dy - np.sin(theta)*Dx)/dist
        cosalpha = (np.cos(theta)*Dx + np.sin(theta)*Dy)/dist

        sinalpha = np.sqrt(np.maximum((1 - cosalpha * cosalpha), 0))
        somme = B1 * sinalpha + B2 * cosalpha*sinalpha
        #-B1 * sinalpha - B2 * cosalpha*sinalpha

        lnEps = np.log(np.maximum(self.lnEps_cr,(dist/a - 2)))
                
        val = (8/5)*self.mu*np.pi*(a**2)*eieijt*somme*lnEps
        val2 = (1/4)*val
        
        return val, val2
        
    def forcesLubrification(self, Dx, Dy, dist, theta):
        #Computes the lubrification forces between two particles
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
        RRi = np.sqrt((xs - self.Nx)**2 + (ys - self.Ny)**2)
        tmp = -6*((self.Es*(self.Nx - xs))/(6*np.pi*self.mu*self.radius*RRi))*(2*(self.radius/RRi)**13-(self.radius/RRi)**7)
        return tmp*xs
    
    def compute_force_squirmer_border_y(self, xs, ys):
        RRi = np.sqrt((xs - self.Nx)**2 + (ys - self.Ny)**2)
        tmp = -6*((self.Es*(self.Ny - ys))/(6*np.pi*self.mu*self.radius*RRi))*(2*(self.radius/RRi)**13-(self.radius/RRi)**7)
        return tmp*ys
    
    def force_torque_lubrification_border_x(self, xs, theta, border):
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
        F_x = -(16/5) * np.pi * self.mu * a * eieijt * somme * lnEps * Dx
        gamma_w = (16/5)*self.mu*np.pi*(a**2)*eieijt*somme*lnEps
        return F_x, gamma_w
    
    def force_torque_lubrification_border_y(self, ys, theta, border):
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
        F_y = -9 * self.mu * np.pi*a*(1/4)*sommeFz* lnEps * Dy
        gamma_w = (16/5)*self.mu*np.pi*(a**2)*eieijt*somme*lnEps
        return F_y, gamma_w

    #Reflective boundary condition
    def ref_border_x(self, xs, orientation, boundary):
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
        if boundary == 1:
            #1 for right boundary
            xs -= 2*self.Nx
        else:
            xs += 2*self.Nx
        return xs

    def loop_time(self):
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


            dist_steric = (dists<self.ds)&(dists!=0)
            dist_lubrification = (dists<=3*a)&(dists!=0)
            j_dist_steric = np.where(dist_steric)
            j_dist_lubr = np.where(dist_lubrification)

            #Steric Forces
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

            #Noise
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
            self.xs += self.dt*(self.v0*np.cos(self.orientations) - self.Fs_x - self.Fs_pw[0] + self.Fl_x) + np.sqrt(2*self.dt*self.Do)*self.nos
            self.ys += self.dt*(self.v0*np.sin(self.orientations) - self.Fs_y - self.Fs_pw[1] + self.Fl_y) + np.sqrt(2*self.dt*self.Do)*self.nos

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
                print(tout)
                data = [self.xs.tolist(), self.ys.tolist(), self.orientations.tolist(),
                        self.Fs_x.tolist(), self.Fs_y.tolist(), self.Fl_x.tolist(), self.Fl_y.tolist(),
                        self.val.tolist(), self.gamma_w.tolist(), self.Fs_pw.tolist(), tout]
                history.append(data)
                tout += self.dt_out

        self.history = history
        return history