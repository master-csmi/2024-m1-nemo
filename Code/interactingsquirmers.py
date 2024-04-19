import matplotlib
import copy
import numpy as np

matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

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
    
    def plot_squirmers_positions(self, history):
        R = self.R
        plt.figure(figsize=(8, 8))
        plt.plot([-R, R], [-R, -R], 'k-', linewidth=2)  # Bottom side
        plt.plot([-R, R], [R, R], 'k-', linewidth=2)  # Top side
        plt.plot([-R, -R], [-R, R], 'k-', linewidth=2)  # Left side
        plt.plot([R, R], [-R, R], 'k-', linewidth=2)  # Right side

        for step in history:
            plt.scatter(step['squirmer1'].x, step['squirmer1'].y, color='blue', s=50)
            #plt.quiver(step['squirmer1'].x, step['squirmer1'].y, np.cos(step['squirmer1'].orientation), np.sin(step['squirmer1'].orientation), color='blue', scale=10)

            plt.scatter(step['squirmer2'].x, step['squirmer2'].y, color='red', s=50)
            #plt.quiver(step['squirmer2'].x, step['squirmer2'].y, np.cos(step['squirmer2'].orientation), np.sin(step['squirmer2'].orientation), color='red', scale=10)
        
        plt.axis('equal')
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.title('Positions and Orientations of Squirmers')
        plt.grid(True)
        plt.show()

    def plot_dist_sq(self, dist_list):
        plt.figure(figsize=(8, 6))
        plt.plot(np.arange(0, self.T-self.dt_out, self.dt_out), dist_list)
        plt.xlabel('Time')
        plt.ylabel('Distance between squirmers')
        plt.title('Distance between squirmers over time')
        plt.grid(True)
        plt.show()

    def plot_dist_border(self, dist_border):
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
        plt.show()

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
        sommeFz = B1*sinalpha*cosalpha - (1/2)*B1*cosalpha*eieijt**2 + B2*sinalpha*cosalpha**2 - (1/2)*B2*(2*cosalpha**2-1)*eieijt**2
        epsilon = max((dist - 2*a)/a, 0.01)

        #lambda = mu = 1
        F_x = np.pi * a * eieijt * somme * np.log(epsilon) * Dx
        F_y = -9*np.pi*a*(1/4)*sommeFz*np.log(epsilon) * Dy

        return F_x, F_y
    
    def loop_time(self):
        tout = self.dt_out
        a = self.squirmer1.radius
        history = []
        dist_list = []
        for t in np.arange(0, self.T, self.dt):
            Fs_x = 0
            Fs_y = 0
            Dx, Dy, dist = self.distance_sq()
            #Force between squirmers
            if dist < self.ds:
                tmp = -3*(self.Es/a)*(Dy/dist)*(2*(2*a/dist)**13-(2*a/dist)**7)
                Fs_x = tmp * Dx
                Fs_y = tmp * Dy
            
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
                sq1_copie = copy.deepcopy(self.squirmer1)
                sq2_copie = copy.deepcopy(self.squirmer2)

                #List that contains positions of squirmers
                history.append({'squirmer1':sq1_copie, 'squirmer2':sq2_copie})

                #List that contains the distance between the squirmers
                dist_list.append(dist)
                tout += self.dt_out

        return history, dist_list
    
    def run(self, dist_sq = False):
        self.init_two_squirmers()
        history, dist_list = self.loop_time()
        self.plot_squirmers_positions(history)
        if (dist_sq == True):
            self.plot_dist_sq(dist_list)