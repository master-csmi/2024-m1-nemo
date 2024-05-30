import numpy as np
from csv_file import export_data_csv, read_csv_file
from plot import plot_dist_sq, plot_squirmers_positions

class InteractingSquirmers:

    def __init__(self, squirmer1, squirmer2, R, dt, dt_out, T, Es, ds, mu, Eo, lnEps_cr):
        self.squirmer1 = squirmer1
        self.squirmer2 = squirmer2
        self.R = R
        self.dt = dt
        self.dt_out = dt_out
        self.Es = Es
        self.mu = mu
        self.Eo = Eo
        self.ds = ds
        self.T = T
        self.lnEps_cr = lnEps_cr
    
    def is_in_square(self):
        #return True if the squirmers are in the square
        return (abs(self.squirmer1.x) < (self.R-self.squirmer1.radius)) and (abs(self.squirmer1.y) < (self.R-self.squirmer1.radius)) and (abs(self.squirmer2.x) < (self.R-self.squirmer2.radius)) and (abs(self.squirmer2.y) < (self.R-self.squirmer2.radius))

    def check_squirmers_square(self):
        #return True if the squirmers have been successfully initialized
        if (self.is_in_square() == False):
            raise ValueError("Squirmers must be inside the square")
        return True
    
    def distance_sq(self):
        #return distance between the two squirmers
        Dx = self.squirmer2.x - self.squirmer1.x
        Dy = self.squirmer2.y - self.squirmer1.y
        return Dx, Dy, np.sqrt(Dx**2 + Dy**2)
    
    def distance_center(self):
        #Compute the distance between the squirmers and the center (0,0)
        dist_sq1 = np.sqrt(self.squirmer1.x**2 + self.squirmer1.y**2)
        dist_sq2 = np.sqrt(self.squirmer2.x**2 + self.squirmer2.y**2)
        return dist_sq1, dist_sq2

    def forcesSteric(self, Dx, Dy, dist):
        a = self.squirmer1.radius
        tmp = -3*(self.Es/a)*(2*(2*a/dist)**13-(2*a/dist)**7)/dist
        Fs_x = tmp * Dx
        Fs_y = tmp * Dy
        return Fs_x, Fs_y
    
    def torquesLubrification(self,choice):
        Dx, Dy, dist = self.distance_sq()
        
        if (choice == 1):
            squirmer = self.squirmer1
        else:
            squirmer = self.squirmer2
            Dx = -Dx
            Dy = -Dy

        theta = squirmer.orientation
        beta = squirmer.beta
        a = squirmer.radius
        
        ex = Dx/dist
        ey = Dy/dist

        lnEps = -np.log(max(self.lnEps_cr,(dist/a - 2)))
                
        val = self.Eo*(1 + beta*(np.cos(theta)*ex + np.sin(theta)*ey))*lnEps*(ex*np.sin(theta) - ey*np.cos(theta))
        
        return val
        
    def forcesLubrification(self, choice):
        Dx, Dy, dist = self.distance_sq()
        if (choice == 1):
            squirmer = self.squirmer1
        else:
            squirmer = self.squirmer2
            Dx = -Dx
            Dy = -Dy

        theta = squirmer.orientation
        B1 = squirmer.B1
        B2 = squirmer.B2
        a = squirmer.radius

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
    
    def compute_force_squirmer_border_x(self, choice):
        if (choice == 1):
            squirmer = self.squirmer1
        else:
            squirmer = self.squirmer2
        x = squirmer.x
        y = squirmer.y
        RRi = np.sqrt((x - self.R)**2 + (y - self.R)**2)
        tmp = 6*((self.Es*(self.R - x))/(squirmer.radius*RRi))*(2*(squirmer.radius/RRi)**13-(squirmer.radius/RRi)**7)
        return tmp*squirmer.x
    
    def compute_force_squirmer_border_y(self, choice):
        if (choice == 1):
            squirmer = self.squirmer1
        else:
            squirmer = self.squirmer2
        y = squirmer.y
        x = squirmer.x
        RRi = np.sqrt((x - self.R)**2 + (y - self.R)**2)
        tmp = 6*((self.Es*(self.R - y))/(squirmer.radius*RRi))*(2*(squirmer.radius/RRi)**13-(squirmer.radius/RRi)**7)
        return tmp*squirmer.y

    def compute_torque_squirmer_border(self, choice, dist_center):
        if (choice == 1):
            squirmer = self.squirmer1
        else:
            squirmer = self.squirmer2
        ex = squirmer.x / dist_center
        ey = squirmer.y / dist_center

        lnEps = -np.log(max(self.lnEps_cr, (self.R - dist_center)/squirmer.radius - 1))

        gamma_w = 2*self.Eo*(1 + squirmer.beta*(np.cos(squirmer.orientation)*ex + np.sin(squirmer.orientation)*ey)) * \
                lnEps*(np.sin(squirmer.orientation)*ex - np.cos(squirmer.orientation)*ey)
        
        return gamma_w
    
    #Reflective boundary condition
    def ref_bound_x(self, choice, boundary):
        if (choice == 1):
            squirmer = self.squirmer1
        else:
            squirmer = self.squirmer2
        x = squirmer.x
        orient_new = np.pi - squirmer.orientation
        diff = self.R - abs(x)
        if boundary == 1:
            #1 for the right border
            x = self.R - diff
        else:
            x = -self.R + diff
        return x, orient_new
    
    def ref_bound_y(self, choice, boundary):
        if (choice == 1):
            squirmer = self.squirmer1
        else:
            squirmer = self.squirmer2
        y = squirmer.y
        orient_new = -squirmer.orientation
        diff = self.R - abs(y)
        if boundary == 1:
            #1 for the up boundary
            y = self.R - diff
        else:
            y = -self.R + diff
        return y, orient_new
    
    def loop_time(self):
        tout = self.dt_out
        a = self.squirmer1.radius
        #List that contains data to export
        history = []
        data = [self.squirmer1.x, self.squirmer1.y, 
                        self.squirmer2.x, self.squirmer2.y, 
                        self.squirmer1.orientation, self.squirmer2.orientation,
                        0, 0, 
                        0, 0, 
                        0, 0]
        history.append(data)

        for t in np.arange(0, self.T, self.dt):
            Fs_x = 0
            Fs_y = 0
            Dx, Dy, dist = self.distance_sq()

            #Force between squirmers
            if dist < self.ds:
                Fs_x, Fs_y = self.forcesSteric(Dx, Dy, dist)
            
            #Compute torques exerted on squirmer by other squirmer
            val1 = 0
            val2 = 0
            
            if (dist < 3*a):
                val1 = self.torquesLubrification(1)
                val2 = self.torquesLubrification(2)
                    
            #Lubrification forces
            Fl_x1 = 0.0
            Fl_y1 = 0.0
            Fl_x2 = 0.0
            Fl_y2 = 0.0
            
            if dist < 3*a:
                Fl_x1, Fl_y1 = self.forcesLubrification(1)
                Fl_x2, Fl_y2 = self.forcesLubrification(2)

            #Force between a squirmer and a border
            Fs_pw1 = [0,0]
            Fs_pw2 = [0,0]
            if ((self.R-abs(self.squirmer1.x)) < 2**(1/6)*self.squirmer1.radius):
                Fs_pw1[0] = self.compute_force_squirmer_border_x(1)
            if ((self.R-abs(self.squirmer1.y)) < 2**(1/6)*self.squirmer1.radius):
                Fs_pw1[1] = self.compute_force_squirmer_border_y(1)
            if ((self.R-abs(self.squirmer2.x)) < 2**(1/6)*self.squirmer1.radius):
                Fs_pw2[0] = self.compute_force_squirmer_border_x(2)
            if ((self.R-abs(self.squirmer2.y)) < 2**(1/6)*self.squirmer1.radius):
                Fs_pw2[1] = self.compute_force_squirmer_border_y(2)
        
            #Evolution of position
            self.squirmer1.orientation += self.dt*(val1 + 0.25*val2)
            self.squirmer1.x += self.dt*(self.squirmer1.velocity * np.cos(self.squirmer1.orientation) + Fs_x + Fl_x1 + Fs_pw1[0])
            self.squirmer1.y += self.dt*(self.squirmer1.velocity * np.sin(self.squirmer1.orientation) + Fs_y + Fl_y1 + Fs_pw1[1])
            
            self.squirmer2.orientation += self.dt*(val2 + 0.25*val1)
            self.squirmer2.x += self.dt*(self.squirmer2.velocity * np.cos(self.squirmer2.orientation) - Fs_x + Fl_x2 + Fs_pw2[0])
            self.squirmer2.y += self.dt*(self.squirmer2.velocity * np.sin(self.squirmer2.orientation) - Fs_y + Fl_y2 + Fs_pw2[1])

            #Compute torque exerted on squirmer by the wall
            gamma_w1 = 0
            gamma_w2 = 0
            dist_sq1, dist_sq2 = self.distance_center()
            if ((self.R - abs(self.squirmer1.x)) < 2**(1/6) * self.squirmer1.radius):
                gamma_w1 += self.compute_torque_squirmer_border(1, dist_sq1)
            if ((self.R - abs(self.squirmer1.y)) < 2**(1/6) * self.squirmer1.radius):
                gamma_w1 += self.compute_torque_squirmer_border(1, dist_sq1)
            if ((self.R - abs(self.squirmer2.x)) < 2**(1/6) * self.squirmer2.radius):
                gamma_w2 += self.compute_torque_squirmer_border(2, dist_sq2)
            if ((self.R - abs(self.squirmer2.y)) < 2**(1/6) * self.squirmer2.radius):
                gamma_w2 += self.compute_torque_squirmer_border(2, dist_sq2)
            
            # if gamma_w1 != 0:
            #     print("gamma_w1 =", gamma_w1)
            # if gamma_w2 != 0:
            #     print("gamma_w2 =", gamma_w2)

            #Update orientation
            self.squirmer2.orientation += self.dt*(gamma_w1)
            self.squirmer2.orientation += self.dt*(gamma_w2)

            #Reflective Boundary
            if ((self.R-self.squirmer1.x) < 2**(1/6)*self.squirmer1.radius):
                self.squirmer1.x, self.squirmer1.orientation = self.ref_bound_x(1,1)
            if ((self.R+self.squirmer1.x) < 2**(1/6)*self.squirmer1.radius):
                self.squirmer1.x, self.squirmer1.orientation = self.ref_bound_x(1,2)
            if ((self.R-self.squirmer1.y) < 2**(1/6)*self.squirmer1.radius):
                self.squirmer1.y, self.squirmer1.orientation = self.ref_bound_y(1,1)
            if ((self.R+self.squirmer1.y) < 2**(1/6)*self.squirmer1.radius):
                self.squirmer1.y, self.squirmer1.orientation = self.ref_bound_y(1,2)

            if ((self.R-self.squirmer2.x) < 2**(1/6)*self.squirmer2.radius):
                self.squirmer2.x, self.squirmer2.orientation = self.ref_bound_x(2,1)
            if ((self.R+self.squirmer2.x) < 2**(1/6)*self.squirmer2.radius):
                self.squirmer2.x, self.squirmer2.orientation = self.ref_bound_x(2,2)
            if ((self.R-self.squirmer2.y) < 2**(1/6)*self.squirmer2.radius):
                self.squirmer2.y, self.squirmer2.orientation = self.ref_bound_y(2,1)
            if ((self.R+self.squirmer2.y) < 2**(1/6)*self.squirmer1.radius):
                self.squirmer2.y, self.squirmer2.orientation = self.ref_bound_y(2,2)
            
            #Update the data to export
            if t >= tout:
                data = [self.squirmer1.x, self.squirmer1.y, 
                        self.squirmer2.x, self.squirmer2.y, 
                        self.squirmer1.orientation, self.squirmer2.orientation,
                        Fl_x1, Fl_y1, 
                        Fl_x2, Fl_y2, 
                        dist, tout]
                history.append(data)
                tout += self.dt_out

        return history
    
    def run(self, file_name_csv, filename_pos='position_graph', filename_dist='dist_squirmer_graph',):
        self.check_squirmers_square()
        history = self.loop_time()
        export_data_csv(file_name_csv, history)
        plot_squirmers_positions(self.R, history, filename_pos)
        plot_dist_sq(history, filename_dist)
