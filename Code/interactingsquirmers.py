import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from csv_file import export_data_csv, read_csv_file
from squirmer import Squirmer
from plot import plot_dist_sq, plot_squirmers_positions

class InteractingSquirmers:

    def __init__(self, N, xs, ys, orientations, radius, beta, v0, R, dt, dt_out, T, Es, ds, mu, Eo, lnEps_cr):
        self.N = N
        self.xs = np.array(xs, dtype=float)
        self.ys = np.array(ys, dtype=float)
        self.orientations = np.array(orientations, dtype=float)
        self.radius = radius
        self.beta = beta
        self.v0 = v0
        self.R = R
        self.dt = dt
        self.dt_out = dt_out
        self.Es = Es
        self.mu = mu
        self.Eo = Eo
        self.ds = ds
        self.T = T
        self.lnEps_cr = lnEps_cr
        self.Fs_x, self.Fs_y, self.Fl_x, self.Fl_y, self.val, self.gamma_w = np.zeros(N, dtype=float), np.zeros(N, dtype=float), np.zeros(N, dtype=float), np.zeros(N, dtype=float), np.zeros(N, dtype=float), np.zeros(N, dtype=float)
        self.Fs_pw = np.zeros((2,N), dtype=float)

        colors = list(matplotlib.colors.CSS4_COLORS.keys())
        colors = [color for color in colors if not self.is_light_color(matplotlib.colors.CSS4_COLORS[color])]
        np.random.shuffle(colors)
        self.colors = colors[:N]
        
        self.squirmers = []
        for i in range(N):
            self.squirmers.append(Squirmer(self.xs[i], self.ys[i], self.orientations[i], radius, beta, v0))

    def is_light_color(self, hex_color):
        #Define what a color too bright is
        rgb = matplotlib.colors.hex2color(hex_color)
        luminance = 0.2126 * rgb[0] + 0.7152 * rgb[1] + 0.0722 * rgb[2]
        return luminance > 0.7
    
    def is_in_square(self):
        #return True if the squirmers are in the square
        return np.all(abs(self.xs) <= (self.R - self.radius)) & np.all(abs(self.ys) <= (self.R - self.radius))

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
    
    def distance_all(self, squirmer):
        Dxs, Dys, dists = [], [], []
        for s in self.squirmers:
            Dx, Dy, dist = self.distance_sq(squirmer, s)
            Dxs.append(Dx)
            Dys.append(Dy)
            dists.append(dist)
        return Dxs, Dys, dists
    
    def distance_center(self, squirmer):
        #Compute the distance between the squirmers and the center (0,0)
        dist = np.sqrt(squirmer.x**2 + squirmer.y**2)
        return dist

    def forcesSteric(self, squirmer1, squirmer2):
        #Compute the steric forces between two particles
        a = self.radius
        Dx, Dy, dist = self.distance_sq(squirmer1, squirmer2)

        tmp = -3*(self.Es/a)*(2*(2*a/dist)**13-(2*a/dist)**7)/np.sqrt(dist)
        Fs_x =  tmp*Dx
        Fs_y = tmp*Dy
        return Fs_x, Fs_y
    
    def torquesLubrification(self, squirmer1, squirmer2):
        #Computes the lubrification torques produced by two interacting squirmers
        Dx, Dy, dist = self.distance_sq(squirmer1, squirmer2)
        
        theta = squirmer1.orientation
        beta = squirmer1.beta
        a = self.radius
        
        ex = Dx/dist
        ey = Dy/dist

        lnEps = -np.log(max(self.lnEps_cr,(dist/a - 2)))
                
        val = self.Eo*(1 + beta*(np.cos(theta)*ex + np.sin(theta)*ey))*lnEps*(ex*np.sin(theta) - ey*np.cos(theta))
        
        return val
        
    def forcesLubrification(self, squirmer1, squirmer2):
        #Computes the lubrification forces between two particles
        Dx, Dy, dist = self.distance_sq(squirmer1, squirmer2)

        theta = squirmer1.orientation
        B1 = squirmer1.B1
        B2 = squirmer1.B2
        a = self.radius

        eieijt = (np.cos(theta)*Dy - np.sin(theta)*Dx)/dist
        cosalpha = (np.cos(theta)*Dx + np.sin(theta)*Dy)/dist

        sinalpha = np.sqrt(1 - cosalpha * cosalpha)
        somme = - B1 * sinalpha - B2 * cosalpha*sinalpha
        sommeFz = B1 * sinalpha * cosalpha - (1/2)*B1 * cosalpha * eieijt**2 + B2 * sinalpha * cosalpha**2 - (1/2)*B2 * (2*cosalpha**2-1) * eieijt**2

        lnEps = -np.log(max(self.lnEps_cr,(dist/a - 2)))
        
        #lambda=1
        F_x = np.pi * self.mu * a * eieijt * somme * lnEps * Dx
        F_y = -9* self.mu * np.pi*a*(1/4)*sommeFz* lnEps * Dy

        return F_x, F_y
    
    def compute_force_squirmer_border_x(self, squirmer):
        x = squirmer.x
        y = squirmer.y
        RRi = np.sqrt((x - self.R)**2 + (y - self.R)**2)
        tmp = -6*((self.Es*(self.R - x))/(squirmer.radius*RRi))*(2*(squirmer.radius/RRi)**13-(squirmer.radius/RRi)**7)
        return tmp*squirmer.x
    
    def compute_force_squirmer_border_y(self, squirmer):
        y = squirmer.y
        x = squirmer.x
        RRi = np.sqrt((x - self.R)**2 + (y - self.R)**2)
        tmp = -6*((self.Es*(self.R - y))/(squirmer.radius*RRi))*(2*(squirmer.radius/RRi)**13-(squirmer.radius/RRi)**7)
        return tmp*squirmer.y

    def compute_torque_squirmer_border(self, squirmer):
        dist_center = np.sqrt(squirmer.x**2 + squirmer.y**2)
        ex = squirmer.x / dist_center
        ey = squirmer.y / dist_center

        lnEps = -np.log(max(self.lnEps_cr, (self.R - dist_center)/squirmer.radius - 1))

        gamma_w = 2*self.Eo*(1 + squirmer.beta*(np.cos(squirmer.orientation)*ex + np.sin(squirmer.orientation)*ey)) * \
                lnEps*(np.sin(squirmer.orientation)*ex - np.cos(squirmer.orientation)*ey)
        
        return gamma_w
    
    #Reflective boundary condition
    def ref_border_x(self, squirmer, boundary):
        diff = abs(self.R - abs(squirmer.x))
        squirmer.orientation = np.pi - squirmer.orientation
        #Keeps orientation between [0, 2pi]
        squirmer.orientation = squirmer.orientation % (2 * np.pi)
        if boundary == 1:
            #1 for the right border
            squirmer.x = self.R - diff
        else:
            squirmer.x = -self.R + diff

        return squirmer.x, squirmer.orientation
    
    def ref_border_y(self, squirmer, boundary):
        #reflective border y
        squirmer.orientation = -squirmer.orientation
        #Keeps orientation between [0, 2pi]
        squirmer.orientation = squirmer.orientation % (2 * np.pi)
        diff = abs(self.R - abs(squirmer.y))
        if boundary == 1:
            #1 for the up boundary
            squirmer.y = self.R - diff
        else:
            squirmer.y = -self.R + diff

        return squirmer.y, squirmer.orientation
    
    def loop_time(self):
        tout = self.dt_out
        a = self.radius
        #List that contains data to export
        history = []
        data = [self.xs.tolist(), self.ys.tolist(), self.orientations.tolist(),
                self.Fs_x.tolist(), self.Fs_y.tolist(), self.Fl_x.tolist(), self.Fl_y.tolist(),
                self.val.tolist(), self.gamma_w.tolist(), self.Fs_pw.tolist()]
        history.append(data)

        for t in np.arange(0, self.T, self.dt):
            # Reset forces and torques
            self.Fs_x.fill(0)
            self.Fs_y.fill(0)
            self.Fl_x.fill(0)
            self.Fl_y.fill(0)
            self.val.fill(0)
            self.gamma_w.fill(0)
            self.Fs_pw.fill(0)
            for i,s in enumerate(self.squirmers):
                Dx, Dy, dist = self.distance_all(s)

                for j in range(len(dist)):
                    #Force between squirmers
                    if (dist[j] != 0) and (dist[j] < self.ds):
                        Fs_x, Fs_y = self.forcesSteric(s, self.squirmers[j])
                        self.Fs_x[i] += Fs_x
                        self.Fs_y[i] += Fs_y
                    
                    #Lubrification forces and torques
                    if (dist[j] != 0) and (dist[j] <= 3*a):
                        #Forces
                        Fl_x, Fl_y = self.forcesLubrification(s, self.squirmers[j])
                        self.Fl_x[i] += Fl_x
                        self.Fl_y[i] += Fl_y
                        self.Fl_x[j] -= Fl_x
                        self.Fl_y[j] -= Fl_y

                        #Torques
                        val1 = self.torquesLubrification(s, self.squirmers[j])
                        val2 = self.torquesLubrification(self.squirmers[j], s)
                        self.val[i] += val1 + 0.25*val2
                        self.val[j] += val2 + 0.25*val1

                #Force between a squirmer and a border
                if ((self.R-abs(s.x)) < 2**(1/6)*a):
                    self.Fs_pw[0][i] += self.compute_force_squirmer_border_x(s)
                if ((self.R-abs(s.y)) < 2**(1/6)*a):
                    self.Fs_pw[1][i] += self.compute_force_squirmer_border_y(s)

                #Torque exerted on squirmer by the wall
                if ((self.R - abs(s.x)) < 2**(1/6) * a):
                    self.gamma_w[i] += self.compute_torque_squirmer_border(s)
                if ((self.R - abs(s.y)) < 2**(1/6) * a):
                    self.gamma_w[i] += self.compute_torque_squirmer_border(s)
        
            #Evolution of position
            self.orientations += self.dt*(self.val + self.gamma_w)
            self.xs += self.dt*(self.v0 * np.cos(self.orientations) - self.Fs_x - self.Fs_pw[0] + self.Fl_x)
            self.ys += self.dt*(self.v0 * np.sin(self.orientations) - self.Fs_y - self.Fs_pw[1] + self.Fl_y)
            for i, s in enumerate(self.squirmers):
                s.x = self.xs[i]
                s.y = self.ys[i]
                s.orientation = self.orientations[i]

                #Reflective Boundary
                if (self.R - s.x) < a:
                    s.x, s.orientation = self.ref_border_x(s, 1)
                    self.xs[i], self.orientations[i] = s.x, s.orientation
                if (self.R + s.x) < a:
                    s.x, s.orientation = self.ref_border_x(s, 2)
                    self.xs[i], self.orientations[i] = s.x, s.orientation
                if (self.R - s.y) < a:
                    s.y, s.orientation = self.ref_border_y(s, 1)
                    self.ys[i], self.orientations[i] = s.y, s.orientation
                if (self.R + s.y) < a:
                    s.y, s.orientation = self.ref_border_y(s, 2)
                    self.ys[i], self.orientations[i] = s.y, s.orientation
            
            #Update the data to export
            if t >= tout:
                data = [self.xs.tolist(), self.ys.tolist(), self.orientations.tolist(),
                        self.Fs_x.tolist(), self.Fs_y.tolist(), self.Fl_x.tolist(), self.Fl_y.tolist(),
                        self.val.tolist(), self.gamma_w.tolist(), self.Fs_pw.tolist()]
                history.append(data)
                tout += self.dt_out

        return history
    
    def run(self, file_name_csv, filename_pos='position_graph', filename_dist='dist_squirmer_graph',):
        self.check_squirmers_square()
        history = self.loop_time()
        export_data_csv(file_name_csv, history)
        plot_squirmers_positions(self.R, history, filename_pos)
        plot_dist_sq(history, filename_dist)

# #velocity
# v0 = 1
# #length of the square
# L = 2
# #half of the length of the square
# R = L/2
# #squirmers' radius
# a = 0.15
# #coordinates
# x1, y1 = -2*a, 0
# x2, y2 = a/1.1, 0
# #orientations
# orient1 = np.pi/2
# orient2 = [-np.pi/4, np.pi/4, np.pi/2, -np.pi/2, 3*np.pi/4, -3*np.pi/4, np.pi, 2*np.pi]
# #border simulation
# xs = [-0.6, 0.6, 0.2, -0.2]
# ys = [-0.2, -0.2, 0.6, -0.6]
# orients = [3*np.pi/4, np.pi/4, 3*np.pi/4, -np.pi/4]
# #files labels
# file = ["-pi.4", "pi.4", "pi.2", "-pi.2", "3pi.4", "-3pi.4", "pi", "2pi"]
# #betas
# betas = [0, -7.5, 7.5]
# beta0, betainf, betasup = 0, -7.5, 7.5
# #time-step
# dt = 1e-4
# #cut-off for -log
# lnEps_cr = np.exp(-5)
# #amplitude of steric interactions
# Es = 1
# #simulation time
# T = 1
# #periodicity of outputs
# dt_out = 0.05
# #amplitude of orientational interactions
# Eo = ((3./10.)*v0/a)
# #distance of steric interactions
# ds = 2**(7./6)*a
# #viscovity parameter
# mu = 0.01

# N = 2
# xs2 = [x1, x2]
# ys2 = [y1, y2]
# orients2 = [orient1, orient2[2]]

# inter = InteractingSquirmers(N, xs2, ys2, orients2, a, beta0, v0, R, dt, dt_out, T, Es, ds, mu, Eo, lnEps_cr)
# inter.loop_time()