import matplotlib
import copy
import numpy as np

matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from squirmer import Squirmer

v0 = 1
L = 2
R = L/2
x1, y1 = -0.2, 0.7
x2, y2 = 0.2, 0.65
orient1, orient2 = 2*np.pi, np.pi
a = 0.05
beta = -7.5
dt = 1e-4
lnEps_cr = 5
Es = 1
T = 1
dt_out = 0.1
Eo = (3/10)*v0/a
ds = 2**(7/6)*a
D = 0
squirmer1 = Squirmer(x1,y1,orient1,a,beta,v0)
squirmer2 = Squirmer(x2,y2,orient2,a,beta,v0)

class InteractingSquirmers:
    
    def __init__(self, squirmer1, squirmer2, R=R, dt=dt, dt_out=dt_out, T=T, Es=Es, ds=ds, Eo=Eo):
        self.squirmer1 = squirmer1
        self.squirmer2 = squirmer2
        self.R = R
        self.dt = dt
        self.dt_out = dt_out
        self.Es = Es
        self.Eo = Eo
        self.ds = ds
        self.T = T

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

    def forcesHydro(self):
        Dx, Dy, dist = self.distance_sq()
        theta = self.squirmer1.orientation
        B1 = self.squirmer1.B1
        B2 = self.squirmer1.B2
        a = self.squirmer1.raidus

        eex = (np.cos(theta)*Dy - np.sin(theta)*Dx)
        eez = (np.cos(theta)*Dx + np.sin(theta)*Dy)/dist
        W1 = -eez
        W2 = (1/2)*(3*(-eez)**2 - 1)
        somme = B1 * W1 + B2 * W2
        sommeFz = B1*W1*eez + (1/2)*B1*eex**2 + (1/2)*B2*(3*W1**2-1)*eez + (3/2)*B2*W1*eex**2
        epsilon = (dist - 2*a)/a

        #lambda = mu = 1
        F_x = -np.pi * a * eex * somme * np.log(epsilon)
        F_z = -9*np.pi*a*(1/4)*sommeFz*np.log(epsilon)
        T_y1 = 0.6 * (a**2) * np.pi * eex * somme * np.log(epsilon)
        T_y2 = 0.4 * np.pi * (a**2) * somme * np.log(epsilon)
        return F_x, F_z, T_y1, T_y2
    
    def compute_force_squirmer_border_x(self, choice):
        if (choice == 1):
            squirmer = self.squirmer1
        else:
            squirmer = self.squirmer2
        x = squirmer.x
        y = squirmer.y
        RRi = np.sqrt((x - R)**2 + (y - R)**2)
        tmp = 6*((self.Es*(self.R - x))/(squirmer.radius*RRi))*(2*(squirmer.radius/RRi)**13-(squirmer.radius/RRi)**7)
        return tmp*squirmer.x
    
    def compute_force_squirmer_border_y(self, choice):
        if (choice == 1):
            squirmer = self.squirmer1
        else:
            squirmer = self.squirmer2
        y = squirmer.y
        x = squirmer.x
        RRi = np.sqrt((x - R)**2 + (y - R)**2)
        tmp = 6*((self.Es*(self.R - y))/(squirmer.radius*RRi))*(2*(squirmer.radius/RRi)**13-(squirmer.radius/RRi)**7)
        return tmp*squirmer.y
        

    #Reflective boundary condition
    def ref_bound_x(self, choice):
        if (choice == 1):
            squirmer = self.squirmer1
        else:
            squirmer = self.squirmer2
        x = squirmer.x
        phi_part = np.arctan2(squirmer.y, squirmer.x)
        x = (self.R - squirmer.radius) * np.cos(phi_part)
        return x
    
    def ref_bound_y(self, choice):
        if (choice == 1):
            squirmer = self.squirmer1
        else:
            squirmer = self.squirmer2
        y = squirmer.y
        phi_part = np.arctan2(squirmer.y, squirmer.x)
        y = (2*(self.R - squirmer.radius) - y) * np.cos(phi_part)
        return y
    
    def loop_time(self):
        tout = dt_out
        history = []
        dist_border = []
        dist_list = []
        print("critère_x =", self.R - self.squirmer1.radius)
        for t in np.arange(0, self.T, self.dt):
            Fs_x = 0
            Fs_y = 0
            Dx, Dy, dist = self.distance_sq()
            #Force between squirmers
            if dist < self.ds:
                tmp = -3*(self.Es/self.squirmer1.radius)*(Dy/dist)*(2*(2*self.squirmer1.radius/dist)**13-(2*self.squirmer1.radius/dist)**7)
                Fs_x = tmp * Dx
                Fs_y = tmp * Dy

            #Force between a squirmer and a border
            Fs_pw1 = [0,0]
            Fs_pw2 = [0,0]
            if ((self.R-self.squirmer1.x) < 2**(1/6)*self.squirmer1.radius):
                Fs_pw1[0] = self.compute_force_squirmer_border_x(1)
            if ((self.R-self.squirmer1.y) < 2**(1/6)*self.squirmer1.radius):
                Fs_pw1[1] = self.compute_force_squirmer_border_y(1)
            if ((self.R-self.squirmer2.x) < 2**(1/6)*self.squirmer1.radius):
                Fs_pw2[0] = self.compute_force_squirmer_border_x(2)
            if ((self.R-self.squirmer2.y) < 2**(1/6)*self.squirmer1.radius):
                Fs_pw2[1] = self.compute_force_squirmer_border_y(2)
            
            #Compute torques exerted on squirmer by other squirmer
            val1 = 0
            val2 = 0
            if dist < 3*a:
                ex = Dx/dist
                ey = Dy/dist

                lnEps = min(lnEps_cr, -np.log(dist/a - 2))
                val1 = Eo * (1 + beta * (np.cos(self.squirmer1.orientation) * ex + np.sin(self.squirmer1.orientation) * ey)) * \
                        lnEps * (ex * np.sin(self.squirmer1.orientation) - ey * np.cos(self.squirmer1.orientation))
                val2 = Eo * (1 + beta * (np.cos(self.squirmer2.orientation) * ex + np.sin(self.squirmer2.orientation) * ey)) * \
                        lnEps * (ex * np.sin(self.squirmer2.orientation) - ey * np.cos(self.squirmer2.orientation))

            #Compute torque exerted on squirmer by the wall
            gamma_w1 = 0
            gamma_w2 = 0
            # if ((R - dist_sq1) < 2**(1/6) * self.squirmer1.radius):
            #     gamma_w1 = compute_torque_squirmer_border(squirmer1, dist_sq1, R, a, Eo, lnEps_cr)
            # if ((R - dist_sq2) < 2**(1/6) * self.squirmer2.radius):
            #     gamma_w2 = compute_torque_squirmer_border(squirmer2, dist_sq2, R, a, Eo, lnEps_cr)
        
            #Evolution of position
            self.squirmer1.x += self.dt*(self.squirmer1.velocity * np.cos(self.squirmer1.orientation) + Fs_x + Fs_pw1[0])
            self.squirmer1.y += self.dt*(self.squirmer1.velocity * np.sin(self.squirmer1.orientation) + Fs_y + Fs_pw1[1])
            #self.squirmer1.orientation += self.dt*(val1 + 0.25*val1 + gamma_w1)

            self.squirmer2.x += self.dt*(self.squirmer2.velocity * np.cos(self.squirmer2.orientation) - Fs_x + Fs_pw2[0])
            self.squirmer2.y += self.dt*(self.squirmer2.velocity * np.sin(self.squirmer2.orientation) - Fs_y + Fs_pw2[1])
            #self.squirmer2.orientation += self.dt*(val2 + 0.25*val2 + gamma_w2)

            #Reflective boundary
            print("squirmer1.x avant =", self.squirmer1.x)
            if abs(self.squirmer1.x) > self.R-self.squirmer1.radius:
                self.squirmer1.x = self.ref_bound_x(1)
            print("squirmer1.x après =", self.squirmer1.x)
            print("squirmer1.y avant =", self.squirmer1.y)
            if abs(self.squirmer1.y) > self.R-self.squirmer1.radius:
                self.squirmer1.y = self.ref_bound_y(1)
            print("squirmer1.y après =", self.squirmer1.y)
            
            print("squirmer2.x avant =", self.squirmer2.x)
            if abs(self.squirmer2.x) > self.R-self.squirmer2.radius:
                self.squirmer2.x = self.ref_bound_x(2)
            print("squirmer2.x après =", self.squirmer2.x)
            print("squirmer2.y avant =", self.squirmer2.y)
            if abs(self.squirmer2.y) > self.R-self.squirmer2.radius:
                self.squirmer2.y = self.ref_bound_y(2)
            print("squirmer2.y après =", self.squirmer2.y, "\n")

            #Plots
            if t >= tout:
                dist_sq1, dist_sq2 = self.distance_center()
                sq1_copie = copy.deepcopy(self.squirmer1)
                sq2_copie = copy.deepcopy(self.squirmer2)
                #List that contains positions of squirmers
                history.append({'squirmer1':sq1_copie, 'squirmer2':sq2_copie})

                #List that contains the distance between the border and each squirmer
                dist_border.append([self.R - dist_sq1, self.R - dist_sq2])

                #List that contains the distance between the squirmers
                dist_list.append(dist)
                tout += dt_out

        return history, dist_list, dist_border
    
    def run(self, dist_sq = False, dist_bord = False):
        self.init_two_squirmers()
        history, dist_list, dist_border = self.loop_time()
        self.plot_squirmers_positions(history)
        if (dist_sq == True):
            self.plot_dist_sq(dist_list)
        if (dist_bord == True):
            self.plot_dist_border(dist_border)

interact_sq = InteractingSquirmers(squirmer1, squirmer2)
history, dist_list, dist_border = interact_sq.loop_time()
interact_sq.run(True, True)