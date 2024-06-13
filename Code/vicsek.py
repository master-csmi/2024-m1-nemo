import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import time
from squirmer import Squirmer

class Vicsek_continous:

    def __init__(self, N, R, L, v0, beta, radius, T, dt, noise):
        self.N = N
        self.R = R
        self.L = L
        self.v0 = v0
        self.T = T
        self.dt = dt
        self.noise = noise
        self.size = L/R
        self.density = (N*R**2)/(L**2)
        self.ratio = v0/R

        orientations = []
        xs = []
        ys = []

        while len(orientations) < N:
            orientation = np.random.uniform(0, 2*np.pi)
            orientations.append(orientation)
        while len(xs) < N:
            x = np.random.uniform(-(L-2*radius)/2, (L-2*radius)/2)
            y = np.random.uniform(-(L-2*radius)/2, (L-2*radius)/2)
            if not any(np.isclose(x, xs, atol=1.5*radius)) and not any(np.isclose(y, ys, atol=1.5*radius)):
                xs.append(x)
                ys.append(y)

        self.particles = []
        for i in range(len(xs)):
            self.particles.append(Squirmer(xs[i], ys[i], orientations[i], radius, beta, v0))

    def ref_border_x(self, particle, boundary):
        #reflective border x
        diff = self.L/2 - abs(particle.x)
        particle.orientation = np.pi - particle.orientation
        if boundary == 1:
            #1 for the right border
            particle.x = self.L/2 - diff
        else:
            particle.x = self.L/2 + diff

        return particle.x, particle.orientation
    
    def ref_border_y(self, particle, boundary):
        #reflective border y
        particle.orientation = -particle.orientation
        diff = self.L/2 - abs(particle.y)
        if boundary == 1:
            #1 for the up boundary
            particle.y = self.L/2 - diff
        else:
            particle.y = -self.L/2 + diff
        return particle.x, particle.orientation

    def distance(self, p1, p2):
        #Compute the distance between two particles
        return np.sqrt((p1.x-p2.x)**2 + (p1.y-p2.y)**2)
    
    def compute_forces(self):
        return 0
    
    def average_orient(self, particle):
        #Compute the average particle's orientation near each particle
        particles_close = [p for p in self.particles if self.distance(p, particle) <= self.R]
        avrg_orient = np.mean([p.orientation for p in particles_close])
        return avrg_orient
    
    def update_orient(self):
        #Update orientation of each particle
        for particle in self.particles:
            avrg_orient = self.average_orient(particle)
            noise = np.random.uniform(-self.noise/2, self.noise/2)
            new_orient = np.arctan2(np.sin(avrg_orient), np.cos(avrg_orient)) + noise
            particle.orientation = new_orient

    def update_position(self):
        #Update position of each particle
        for particle in self.particles:
            particle.x += self.v0*self.dt*np.cos(particle.orientation)
            particle.y += self.v0*self.dt*np.sin(particle.orientation)
            print(particle.x)

            if particle.x >= self.L/2:
                particle.x, particle.orientation = self.ref_border_x(particle, 1)
            if particle.x <= -self.L/2:
                particle.x, particle.orientation = self.ref_border_x(particle, 2)
            if particle.y >= self.L/2:
                particle.y, particle.orientation = self.ref_border_y(particle, 1)
            if particle.y <= -self.L/2:
                particle.y, particle.orientation = self.ref_border_y(particle, 2)

    def loop_time(self):
        for _ in np.arange(self.dt, self.T, self.dt):
            self.update_orient()
            self.update_position()

    def plot(self, ax):
        ax.scatter([p.x for p in self.particles], [p.y for p in self.particles])
        ax.set_xlim(-self.L / 2, self.L / 2)
        ax.set_ylim(-self.L / 2, self.L / 2)
        ax.set_aspect('equal')

# Exemple d'utilisation
N = 10
R = 1.0
L = 10.0
v0 = 1.0
beta = 0.5
radius = 0.1
T = 10.0
dt = 1
noise = 0.1

vicsek_model = Vicsek_continous(N, R, L, v0, beta, radius, T, dt, noise)

# Plot initial positions and save the figure
fig, ax = plt.subplots(figsize=(8, 8))
vicsek_model.plot(ax)
plt.title("Initial Positions")
plt.savefig("vicsek_initial_positions.png")
plt.close()

# Run the simulation and plot at intervals
num_steps = 2
for step in range(num_steps):
    start_time = time.time()
    vicsek_model.loop_time()
    end_time = time.time()
    sim_time = end_time - start_time
    print(f"Simulation {step + 1} took {sim_time:.2f} seconds")
    fig, ax = plt.subplots(figsize=(8, 8))
    vicsek_model.plot(ax)
    plt.title(f"Positions at Step {step + 1}")
    plt.savefig(f"vicsek_positions_step_{step + 1}.png")
    plt.close()