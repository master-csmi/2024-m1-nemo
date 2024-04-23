import matplotlib

import copy
import os
import numpy as np
import os
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from csv_file import export_data_csv, read_csv_file

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
    
    def plot_squirmers_positions(self, history, filename, dir='graphs'):
        #Plot the position of each squirmers
        R = self.R
        plt.figure(figsize=(8, 8))
        plt.plot([-R, R], [-R, -R], 'k-', linewidth=2)
        plt.plot([-R, R], [R, R], 'k-', linewidth=2)
        plt.plot([-R, -R], [-R, R], 'k-', linewidth=2)
        plt.plot([R, R], [-R, R], 'k-', linewidth=2)

        squirmer1_x = []
        squirmer1_y = []
        squirmer1_orient = []
        squirmer2_x = []
        squirmer2_y = []
        squirmer2_orient = []
        time = []

        for step in history:
            squirmer1_x.append(step[0])
            squirmer1_y.append(step[1])
            squirmer1_orient.append(step[4])
            squirmer2_x.append(step[2])
            squirmer2_y.append(step[3])
            squirmer2_orient.append(step[5])
            time.append(step[-1])
        
        #Squirmers
        plt.scatter(squirmer1_x, squirmer1_y, color='blue', s=10, label = 'Squirmer 1')
        plt.scatter(squirmer2_x, squirmer2_y, color='red', s=10, label= 'Squirmer 2')
        
        for i in range(len(squirmer2_orient)):
            #Orientation
            plt.quiver(squirmer2_x[i], squirmer2_y[i], np.cos(squirmer2_orient[i]), np.sin(squirmer2_orient[i]), color='red', scale=20, width=0.002)
            plt.quiver(squirmer1_x[i], squirmer1_y[i], np.cos(squirmer1_orient[i]), np.sin(squirmer1_orient[i]), color='blue', scale=20, width=0.002)

        plt.axis('equal')
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.title('Positions and Orientations of Squirmers')
        plt.legend()
        plt.grid(True)
        
        if not os.path.exists(dir):
            os.makedirs(dir)
        save_path = os.path.join(dir, filename + '.png')
        plt.savefig(save_path)

    def plot_dist_sq(self, history, filename, dir='graphs'):
        plt.figure(figsize=(8, 6))
        dist_list = []
        time_list = []
        for step in history:
            dist_list.append(step[10])
            time_list.append(step[11])

        plt.plot(time_list, dist_list, label="Distance")
        
        plt.xlabel('Time')
        plt.ylabel('Distance between squirmers')
        plt.title('Distance between squirmers over time')
        plt.grid(True)

        if not os.path.exists(dir):
            os.makedirs(dir)
        save_path = os.path.join(dir, filename + '.png')
        plt.savefig(save_path)

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
                
        val = self.Eo * (1 + beta * (np.cos(theta) * ex + np.sin(theta) * ey)) * lnEps * (ex * np.sin(theta) - ey * np.cos(theta))
        
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
    
    def loop_time(self):
        tout = self.dt_out
        a = self.squirmer1.radius
        #List that contains data to export
        history = []
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
            
            if dist < 3*a:
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

            #Evolution of position
            self.squirmer1.orientation += self.dt*(val1 + 0.25*val2)
            self.squirmer1.x += self.dt*(self.squirmer1.velocity * np.cos(self.squirmer1.orientation) + Fs_x + Fl_x1)
            self.squirmer1.y += self.dt*(self.squirmer1.velocity * np.sin(self.squirmer1.orientation) + Fs_y + Fl_y1)
            
            self.squirmer2.orientation += self.dt*(val2 + 0.25*val1)
            self.squirmer2.x += self.dt*(self.squirmer2.velocity * np.cos(self.squirmer2.orientation) - Fs_x + Fl_x2)
            self.squirmer2.y += self.dt*(self.squirmer2.velocity * np.sin(self.squirmer2.orientation) - Fs_y + Fl_y2)
            
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
        self.plot_squirmers_positions(history, filename_pos)
        self.plot_dist_sq(history, filename_dist)
