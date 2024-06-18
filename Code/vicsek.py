import sys
import os
import copy
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import time
from squirmer import Squirmer

class Vicsek_continous:

    def __init__(self, N, R, L, v0, beta, radius, Es, ds, mu, Eo, lnEps_cr, T, dt, noise):
        self.N = N
        self.R = R
        self.L = L
        self.v0 = v0
        self.T = T
        self.dt = dt
        self.noise = noise
        self.radius = radius
        self.Es = Es
        self.ds = ds
        self.mu = mu
        self.Eo = Eo
        self.lnEps_cr = lnEps_cr
        self.size = L/R
        self.density = (N*R**2)/(L**2)
        self.ratio = v0/R

        self.orientations = np.zeros(N)
        self.xs = np.zeros(N)
        self.ys = np.zeros(N)
        colors = list(matplotlib.colors.CSS4_COLORS.keys())
        colors = [color for color in colors if not self.is_light_color(matplotlib.colors.CSS4_COLORS[color])]
        np.random.shuffle(colors)
        self.colors = colors[:N]

        for j in range(N):
            orientation = np.random.uniform(0, 2*np.pi)
            self.orientations[j] = orientation
        k = 0
        while k < N:
            x = np.random.uniform(-(L-2*radius)/2, (L-2*radius)/2)
            y = np.random.uniform(-(L-2*radius)/2, (L-2*radius)/2)
            #Each particle must have a unique initial position
            if not any(np.isclose(x, self.xs, atol=1.5*radius)) and not any(np.isclose(y, self.ys, atol=1.5*radius)):
                self.xs[k] = x
                self.ys[k] = y
                k += 1

        self.particles = []
        for i in range(len(self.xs)):
            self.particles.append(Squirmer(self.xs[i], self.ys[i], self.orientations[i], radius, beta, v0))

    def dist_particles(self, particle):
        #Compute the distance between the particle in argument and the other ones
        dist = [self.distance(particle, p) for p in self.particles]
        return dist
    
    def how_many_in_square(self):
        #Print and returns the percentage of particle inside the square
        n = 0
        for particle in self.particles:
            if abs(particle.x) <= self.L/2 and abs(particle.y) <=self.L/2:
                n+=1
        percentage = (n/self.N)*100
        print(f"percentage of particles in the square : {percentage}%")
        return percentage

    def is_light_color(self, hex_color):
        #Define what a color too bright is
        rgb = matplotlib.colors.hex2color(hex_color)
        luminance = 0.2126 * rgb[0] + 0.7152 * rgb[1] + 0.0722 * rgb[2]
        return luminance > 0.7

    def polar_order_parameter(self):
        #Returns the polar order parameter
        summ = abs(sum(v0*(np.cos(particle.orientation) + np.sin(particle.orientation)) for particle in self.particles))
        polar_order = 1/(self.N*self.v0)*summ
        return polar_order

    def ref_border_x(self, particle, boundary):
        #reflective border x
        diff = abs(self.L/2 - abs(particle.x))
        particle.orientation = np.pi - particle.orientation
        #Keeps orientation between [0, 2pi]
        particle.orientation = particle.orientation % (2 * np.pi)
        if boundary == 1:
            #1 for the right border
            particle.x = self.L/2 - diff
        else:
            particle.x = -self.L/2 + diff

        return particle.x, particle.orientation
    
    def ref_border_y(self, particle, boundary):
        #reflective border y
        particle.orientation = -particle.orientation
        #Keeps orientation between [0, 2pi]
        particle.orientation = particle.orientation % (2 * np.pi)
        diff = abs(self.L/2 - abs(particle.y))
        if boundary == 1:
            #1 for the up boundary
            particle.y = self.L/2 - diff
        else:
            particle.y = -self.L/2 + diff

        return particle.y, particle.orientation

    def distance(self, p1, p2):
        #Compute the distance between two particles
        return np.sqrt((p1.x-p2.x)**2 + (p1.y-p2.y)**2)
    
    def vector_x(self, p1, p2):
        return p2.x - p1.x
    
    def vector_y(self, p1, p2):
        return p2.y - p1.y
    
    def average_orient(self, particle):
        #Compute the average particle's orientation near each particle
        particles_close = [p for p in self.particles if self.distance(p, particle) <= self.R]
        avrg_orient = np.mean([p.orientation for p in particles_close])
        return avrg_orient
    
    def update_orient(self):
        #Update orientation of each particle
        for i, particle in enumerate(self.particles):
            particle.orientation = self.orientations[i]

        for i, particle in enumerate(self.particles):
            avrg_orient = self.average_orient(particle)
            noise = np.random.uniform(-self.noise/2, self.noise/2)
            new_orient = np.arctan2(np.sin(avrg_orient), np.cos(avrg_orient)) + noise
            particle.orientation = new_orient
            self.orientations[i] = new_orient

    def update_position(self):
        #Update position of each particle
        self.xs += self.v0*self.dt*(np.cos(self.orientations))
        self.ys += self.v0*self.dt*(np.sin(self.orientations))
        for i, particle in enumerate(self.particles):
            particle.x = self.xs[i]
            particle.y = self.ys[i]

            if self.L/2 <= particle.x + self.radius:
                tmpx, tmporientation = self.ref_border_x(particle, 1)
                self.xs[i], particle.x = tmpx, tmpx
                self.orientations[i], particle.orientation = tmporientation, tmporientation
            if -self.L/2 >= particle.x - self.radius:
                tmpx, tmporientation = self.ref_border_x(particle, 2)
                self.xs[i], particle.x = tmpx, tmpx
                self.orientations[i], particle.orientation = tmporientation, tmporientation

            if self.L/2 <= particle.y + self.radius:
                tmpy, tmporientation = self.ref_border_y(particle, 1)
                self.ys[i], particle.y = tmpy, tmpy
                self.orientations[i], particle.orientation = tmporientation, tmporientation
            if -self.L/2 >= particle.y - self.radius:
                tmpy, tmporientation = self.ref_border_y(particle, 2)
                self.ys[i], particle.y = tmpy, tmpy
                self.orientations[i], particle.orientation = tmporientation, tmporientation

    def loop_time(self):
        #Compute the motion of the particles
        for _ in np.arange(self.dt, self.T, self.dt):
            self.update_orient()
            self.update_position()

    def plot(self, ax):
        ax.scatter([p.x for p in self.particles], [p.y for p in self.particles], color=[color for color in self.colors])
        ax.quiver([p.x for p in self.particles], [p.y for p in self.particles], 
                  [np.cos(p.orientation) for p in self.particles], [np.sin(p.orientation) for p in self.particles],
                  color=[color for color in self.colors])
        ax.set_xlim(-self.L / 2, self.L / 2)
        ax.set_ylim(-self.L / 2, self.L / 2)
        ax.set_aspect('equal')

N = 20
R = 0.25
L = 10.0
v0 = 1.0
beta = 0.5
radius = 0.1
T = 1
dt = 0.1
noise = 1e-4
Es = 1
ds = 2**(7./6)*radius
Eo = ((3./10.)*v0/radius)
mu = 0.01
lnEps_cr = np.exp(-5)

vicsek_model = Vicsek_continous(N, R, L, v0, beta, radius, Es, ds, mu, Eo, lnEps_cr, T, dt, noise)

#Plots initial positions and save the figure
fig, ax = plt.subplots(figsize=(8, 8))
vicsek_model.plot(ax)
plt.title("Initial Positions")
plt.savefig("vicsek_initial_positions.png")
plt.close()
# w = 0
# for particle in vicsek_model.particles:
#     print(f"x{w} = {particle.x}")
#     print(f"y{w} = {particle.y}")
#     w+=1
polar = vicsek_model.polar_order_parameter()
print(f"polar parameter = {polar}")
prct = 100
i = 0
compare = 0
#Runs the simulation and plot at intervals
num_steps = 10
for step in range(num_steps):
# while prct == 100:
    compare = copy.deepcopy(vicsek_model.particles)
    start_time = time.time()
    vicsek_model.loop_time()
    end_time = time.time()
    sim_time = end_time - start_time
    print(f"Simulation {i + 1} took {sim_time:.2f} seconds")
    prct = vicsek_model.how_many_in_square()
    i+=1
    # w = 0
    # for particle in vicsek_model.particles:
    #     print(f"x{w} = {particle.x}")
    #     print(f"y{w} = {particle.y}")
    #     w+=1

    polar = vicsek_model.polar_order_parameter()
    print(f"polar parameter = {polar}")

    fig, ax = plt.subplots(figsize=(8, 8))
    vicsek_model.plot(ax)
    plt.title(f"Positions at Step {i + 1}")
    plt.savefig(f"vicsek_positions_step_{i + 1}.png")
    plt.close()