import numpy as np
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from csv_file import export_data_csv, read_csv_file
from squirmer import Squirmer
from plot import plot_squirmers_positions

class InteractingSquirmers:

    def __init__(self, N, xs, ys, orientations, radius, beta, v0, Nx, Ny, dt, dt_out, T, Es, ds, mu, Eo, lnEps_cr, Do, no, border=True):
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
        self.Eo = Eo
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
        Dxs = self.xs - self.xs[:, None]
        Dys = self.ys - self.ys[:, None]
        dists = np.sqrt(Dxs**2 + Dys**2)
        return Dxs, Dys, dists
    
    def distance_center(self, squirmer):
        #Compute the distance between the squirmers and the center (0,0)
        dist = np.sqrt(squirmer.x**2 + squirmer.y**2)
        return dist

    def forcesSteric(self, Dxs, Dys, dists):
        #Compute the steric forces between two particles
        a = self.radius

        tmp = -3*(self.Es/a)*(2*(2*a/dists)**13-(2*a/dists)**7)/np.sqrt(dists)
        Fs_x =  tmp*Dxs
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
        somme = - B1 * sinalpha - B2 * cosalpha*sinalpha

        lnEps = -np.log(np.maximum(self.lnEps_cr,(dist/a - 2)))
                
        val = (16/10)*self.mu*np.pi*(a**2)*eieijt*somme*lnEps
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
        somme = - B1 * sinalpha - B2 * cosalpha*sinalpha
        sommeFz = B1 * sinalpha * cosalpha - (1/2)*B1 * cosalpha * eieijt**2 + B2 * sinalpha * cosalpha**2 - (1/2)*B2 * (2*cosalpha**2-1) * eieijt**2

        lnEps = -np.log(np.maximum(self.lnEps_cr,(dist/a - 2)))
        
        #lambda=1
        F_x = np.pi * self.mu * a * eieijt * somme * lnEps * Dx
        F_y = -9* self.mu * np.pi*a*(1/4)*sommeFz* lnEps * Dy

        return F_x, F_y
    
    def compute_force_squirmer_border_x(self, xs, ys):
        RRi = np.sqrt((xs - self.Nx)**2 + (ys - self.Ny)**2)
        tmp = -6*((self.Es*(self.Nx - xs))/(self.radius*RRi))*(2*(self.radius/RRi)**13-(self.radius/RRi)**7)
        return tmp*xs
    
    def compute_force_squirmer_border_y(self, xs, ys):
        RRi = np.sqrt((xs - self.Nx)**2 + (ys - self.Ny)**2)
        tmp = -6*((self.Es*(self.Ny - ys))/(self.radius*RRi))*(2*(self.radius/RRi)**13-(self.radius/RRi)**7)
        return tmp*ys

    def compute_torque_squirmer_border(self, xs, ys, orientations):
        Dy = self.Ny-abs(ys)
        Dx = self.Nx-abs(xs)
        dist = np.sqrt(Dx**2 + Dy**2)
        theta = orientations
        B1 = self.B1
        B2 = self.B2
        a = self.radius

        eieijt = (np.cos(theta)*Dy - np.sin(theta)*Dx)/dist
        cosalpha = (np.cos(theta)*Dx + np.sin(theta)*Dy)/dist

        sinalpha = np.sqrt(np.maximum((1 - cosalpha * cosalpha), 0))
        somme = - B1 * sinalpha - B2 * cosalpha*sinalpha
        lnEps = -np.log(np.maximum(self.lnEps_cr,(dist/a - 2)))

        gamma_w = (16/5)*self.mu*np.pi*(a**2)*eieijt*somme*lnEps
        # dist_center = np.sqrt(squirmer.x**2 + squirmer.y**2)
        # ex = squirmer.x / dist_center
        # ey = squirmer.y / dist_center

        # lnEps = -np.log(max(self.lnEps_cr, (self.R - dist_center)/squirmer.radius - 1))

        # gamma_w = 2*self.Eo*(1 + squirmer.beta*(np.cos(squirmer.orientation)*ex + np.sin(squirmer.orientation)*ey)) * \
        #         lnEps*(np.sin(squirmer.orientation)*ex - np.cos(squirmer.orientation)*ey)
        
        return gamma_w
    
    #Reflective boundary condition
    def ref_border_x(self, xs, orientation, boundary):
        diff = abs(self.Nx - abs(xs))
        orientation = np.pi - orientation
        #Keeps orientation between [0, 2pi]
        # squirmer.orientation = squirmer.orientation % (2 * np.pi)
        if boundary == 1:
            #1 for the right border
            xs = self.Nx - diff
        else:
            xs = -self.Nx + diff

        return xs, orientation
    
    def ref_border_y(self, ys, orientation, boundary):
        orientation = -orientation
        #Keeps orientation between [0, 2pi]
        # squirmer.orientation = squirmer.orientation % (2 * np.pi)
        diff = abs(self.Ny - abs(ys))
        if boundary == 1:
            #1 for the up boundary
            ys = self.Ny - diff
        else:
            ys = -self.Ny + diff

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
                min_dist = min(dist_nz.flatten() - 2*a)
                self.vector_dists_min.append(min_dist)

            for i, s in enumerate(self.squirmers):
                dist_steric = (dists[i,:]<self.ds)&(dists[i,:]!=0)
                dist_lubrification = (dists[i,:]<=3*a)&(dists[i,:]!=0)

                #Steric forces
                if np.any(dist_steric):
                    j_dist_steric = np.where(dist_steric)[0]
                    Fs_x, Fs_y = self.forcesSteric(Dxs[i, j_dist_steric], Dys[i, j_dist_steric], dists[i, j_dist_steric])
                    self.Fs_x[i] -= np.sum(Fs_x)
                    self.Fs_y[i] -= np.sum(Fs_y)

                #Lubrification forces and torques
                if np.any(dist_lubrification):
                    j_dist_lubr = np.where(dist_lubrification)[0]
    
                    # Calcul vectorisé des forces et torques pour toutes les paires
                    Dx_lubr = Dxs[i, j_dist_lubr]
                    Dy_lubr = Dys[i, j_dist_lubr]
                    dist_lubr = dists[i, j_dist_lubr]
                    orientations_lubr = self.orientations[j_dist_lubr]
                    
                    Fl_x, Fl_y = self.forcesLubrification(Dx_lubr, Dy_lubr, dist_lubr, orientations_lubr)
                    val1, val2 = self.torquesLubrification(Dx_lubr, Dy_lubr, dist_lubr, orientations_lubr)
                    
                    # Mise à jour des forces et des torques
                    self.Fl_x[i] += np.sum(Fl_x)
                    self.Fl_y[i] += np.sum(Fl_y)
                    self.Fl_x[j_dist_lubr] -= Fl_x
                    self.Fl_y[j_dist_lubr] -= Fl_y
                    self.val[i] += np.sum(val1)
                    self.val[j_dist_lubr] += val2

                #Noise
                self.nos[i] = np.random.uniform(-self.no/2, self.no/2)

            #Force between a squirmer and a border
            dist_forces_x = (self.Nx-abs(self.xs)) < 2**(1/6)*a
            dist_forces_y = (self.Ny-abs(self.ys)) < 2**(1/6)*a

            dist_torques_x = (self.Nx - abs(self.xs)) < 2*a
            dist_torques_y = (self.Ny - abs(self.ys)) < 2*a

            #If we simulate in a box
            if self.border:
                i_dist_force_x = np.where(dist_forces_x)[0]
                i_dist_torque_x = np.where(dist_torques_x)[0]
                self.Fs_pw[0][i_dist_force_x] += self.compute_force_squirmer_border_x(self.xs[i_dist_force_x], self.ys[i_dist_force_x])
                self.gamma_w[i_dist_torque_x] += self.compute_torque_squirmer_border(self.squirmers[i_dist_torque_x], self.ys[i_dist_torque_x], self.orientations[i_dist_torque_x])

            i_dist_force_y = np.where(dist_forces_y)[0]
            i_dist_torque_y = np.where(dist_torques_y)[0]
            self.Fs_pw[1][i_dist_force_y] += self.compute_force_squirmer_border_y(self.xs[i_dist_force_y], self.ys[i_dist_force_y])
            self.gamma_w[i_dist_torque_y] += self.compute_torque_squirmer_border(self.xs[i_dist_torque_y], self.ys[i_dist_torque_y], self.orientations[i_dist_torque_y])

            self.orientations += self.dt*(self.val + self.gamma_w) + np.sqrt(2*self.dt*self.Do)*self.nos
            self.xs += self.dt*(self.v0*np.cos(self.orientations) - self.Fs_x - self.Fs_pw[0] + self.Fl_x)
            self.ys += self.dt*(self.v0*np.sin(self.orientations) - self.Fs_y - self.Fs_pw[1] + self.Fl_y)
            
            self.list_polar.append(self.polar_order_parameter())

            #Borders
            mask_x1 = (self.Nx - self.xs) <= a
            mask_x2 = (self.Nx + self.xs) <= a
            mask_y1 = (self.Ny - self.ys) <= a
            mask_y2 = (self.Ny + self.ys) <= a
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
            
            if t >= tout:
                print(tout)
                data = [self.xs.tolist(), self.ys.tolist(), self.orientations.tolist(),
                        self.Fs_x.tolist(), self.Fs_y.tolist(), self.Fl_x.tolist(), self.Fl_y.tolist(),
                        self.val.tolist(), self.gamma_w.tolist(), self.Fs_pw.tolist(), tout]
                history.append(data)
                tout += self.dt_out
                # print(f"polar parameter : {polar}")

        self.history = history
        return history


    def run(self, file_name_csv, filename_pos='position_graph', filename_dist='dist_squirmer_graph',):
        self.check_squirmers_square()
        history = self.loop_time()
        export_data_csv(file_name_csv, history)
        plot_squirmers_positions(self.R, history, filename_pos)