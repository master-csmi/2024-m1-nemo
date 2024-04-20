import matplotlib
import numpy as np
import os
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from csv_file import export_data_csv, read_csv_file

class InteractingSquirmers:

    def __init__(self, squirmer1, squirmer2, R, dt, dt_out, T, Es, ds, Eo, lnEps_cr):
        self.squirmer1 = squirmer1
        self.squirmer2 = squirmer2
        self.R = R
        self.dt = dt
        self.dt_out = dt_out
        self.Es = Es
        self.Eo = Eo
        self.ds = ds
        self.T = T
        self.lnEps_cr = lnEps_cr

    def distance_center(self):
        #Compute the distance between the squirmers and the center (0,0)
        dist_sq1 = np.sqrt(self.squirmer1.x**2 + self.squirmer1.y**2)
        dist_sq2 = np.sqrt(self.squirmer2.x**2 + self.squirmer2.y**2)
        return dist_sq1, dist_sq2
    
    def is_in_square(self):
        #return True if the squirmers are in the square
        return (abs(self.squirmer1.x) < (self.R-self.squirmer1.radius)) and (abs(self.squirmer1.y) < (self.R-self.squirmer1.radius)) and (abs(self.squirmer2.x) < (self.R-self.squirmer2.radius)) and (abs(self.squirmer2.y) < (self.R-self.squirmer2.radius))

    def init_two_squirmers(self):
        #return True if the squirmers have been successfully initialized
        if (self.is_in_square() == False):
            raise ValueError("Squirmers must be inside the square")
        return True
    
    def distance_sq(self):
        #return distance between the two squirmers
        Dx = self.squirmer2.x - self.squirmer1.x
        Dy = self.squirmer2.y - self.squirmer1.y
        return Dx, Dy, np.sqrt(Dx**2 + Dy**2)
    
    def plot_squirmers_positions(self, history, filename='position_graph', dir='graphs'):
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

        for step in history:
            squirmer1_x.append(step[0])
            squirmer1_y.append(step[1])
            squirmer1_orient.append(step[4])
            squirmer2_x.append(step[2])
            squirmer2_y.append(step[3])
            squirmer2_orient.append(step[5])
        
        plt.scatter(squirmer1_x, squirmer1_y, color='blue', s=10, label = 'Squirmer 1')
        plt.scatter(squirmer2_x, squirmer2_y, color='red', s=10, label= 'Squirmer 2')
        
        for i in range(len(squirmer2_orient)):
            plt.quiver(squirmer2_x[i], squirmer2_y[i], np.cos(squirmer2_orient[i]), np.sin(squirmer2_orient[i]), color='red', width=0.002)
            plt.quiver(squirmer1_x[i], squirmer1_y[i], np.cos(squirmer1_orient[i]), np.sin(squirmer1_orient[i]), color='blue', width=0.002)


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

    def plot_dist_sq(self, history, filename='dist_squirmer_graph', dir='graphs'):
        plt.figure(figsize=(8, 6))
        dist_list = []
        for step in history:
            dist_list.append(step[9])
        plt.plot(np.arange(0, self.T-self.dt_out, self.dt_out), dist_list)
        plt.xlabel('Time')
        plt.ylabel('Distance between squirmers')
        plt.title('Distance between squirmers over time')
        plt.grid(True)

        if not os.path.exists(dir):
            os.makedirs(dir)
        save_path = os.path.join(dir, filename + '.png')
        plt.savefig(save_path)

    def plot_dist_border(self, dist_border, filename='dist_to_border_graph', dir='graphs'):
        dist_x = [item[0] for item in dist_border]
        dist_y = [item[1] for item in dist_border]
        plt.figure(figsize=(8, 6))
        plt.plot(np.arange(0, self.T-self.dt_out, self.dt_out), dist_x, label='Distance between squirmer1 and border')
        plt.plot(np.arange(0, self.T-self.dt_out, self.dt_out), dist_y, label='Distance between squirmer2 and border')
        plt.xlabel('Time')
        plt.ylabel('Distance to border')
        plt.title('Distance to border over time')
        plt.legend()
        plt.grid(True)
        
        if not os.path.exists(dir):
            os.makedirs(dir)
        save_path = os.path.join(dir, filename + '.png')
        plt.savefig(save_path)

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

        print(Dx, Dy, dist)

        eieijt = (np.cos(theta)*Dy - np.sin(theta)*Dx)/dist
        cosalpha = (np.cos(theta)*Dx + np.sin(theta)*Dy)/dist
        print(cosalpha, "\n")
        sinalpha = np.sqrt(1 - cosalpha * cosalpha)
        somme = - B1 * sinalpha - B2 * cosalpha*sinalpha
        sommeFz = B1*sinalpha*cosalpha - (1/2)*B1*cosalpha*eieijt**2 + B2*sinalpha*cosalpha**2 - (1/2)*B2*(2*cosalpha**2-1)*eieijt**2
        epsilon = max((dist - 2*a)/a, 0.01)

        #lambda = mu = 1
        F_x = np.pi * a * eieijt * somme * np.log(epsilon) * Dx
        F_y = -9*np.pi*a*(1/4)*sommeFz*np.log(epsilon) * Dy

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
                tmp_y = -3*(self.Es/a)*(Dy/dist)*(2*(2*a/dist)**13-(2*a/dist)**7)
                tmp_x = -3*(self.Es/a)*(Dx/dist)*(2*(2*a/dist)**13-(2*a/dist)**7)
                Fs_x = tmp_x * Dx
                Fs_y = tmp_y * Dy
            
            #Compute torques exerted on squirmer by other squirmer
            val1 = 0
            val2 = 0
            if dist < 3*a:
                ex = Dx/dist
                ey = Dy/dist

                lnEps = min(self.lnEps_cr, -np.log(dist/a - 2))
                val1 = self.Eo * (1 + self.squirmer1.beta * (np.cos(self.squirmer1.orientation) * ex + np.sin(self.squirmer1.orientation) * ey)) * \
                        lnEps * (ex * np.sin(self.squirmer1.orientation) - ey * np.cos(self.squirmer1.orientation))
                val2 = self.Eo * (1 + self.squirmer2.beta * (np.cos(self.squirmer2.orientation) * ex + np.sin(self.squirmer2.orientation) * ey)) * \
                        lnEps * (ex * np.sin(self.squirmer2.orientation) - ey * np.cos(self.squirmer2.orientation))

            #Lubrification forces
            Fl_x1, Fl_y1 = self.forcesLubrification(1)
            Fl_x2, Fl_y2 = self.forcesLubrification(2)

            #Evolution of position
            self.squirmer1.x += self.dt*(self.squirmer1.velocity * np.cos(self.squirmer1.orientation) + Fs_x + Fl_x1)
            self.squirmer1.y += self.dt*(self.squirmer1.velocity * np.sin(self.squirmer1.orientation) + Fs_y + Fl_y1)
            self.squirmer1.orientation += self.dt*(val1 + 0.25*val1)

            self.squirmer2.x += self.dt*(self.squirmer2.velocity * np.cos(self.squirmer2.orientation) - Fs_x + Fl_x2)
            self.squirmer2.y += self.dt*(self.squirmer2.velocity * np.sin(self.squirmer2.orientation) - Fs_y + Fl_y2)
            self.squirmer2.orientation += self.dt*(val2 + 0.25*val2)

            #Plots
            if t >= tout:
                data = [self.squirmer1.x, self.squirmer1.y, 
                        self.squirmer2.x, self.squirmer2.y, 
                        self.squirmer1.orientation, self.squirmer2.orientation,
                        Fl_x1, Fl_y1, 
                        Fl_x2, Fl_y2, 
                        dist]
                history.append(data)
                tout += self.dt_out

        return history
    
    def run(self, file_name):
        self.init_two_squirmers()
        history = self.loop_time()
        export_data_csv(file_name, history)
        self.plot_squirmers_positions(history)
        self.plot_dist_sq(history)
