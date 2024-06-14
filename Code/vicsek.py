import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import time
from squirmer import Squirmer

class Vicsek_continous:

    def __init__(self, N, R, L, v0, beta, radius, Es, T, dt, noise):
        self.N = N
        self.R = R
        self.L = L
        self.v0 = v0
        self.T = T
        self.dt = dt
        self.noise = noise
        self.radius = radius
        self.Es = Es
        self.size = L/R
        self.density = (N*R**2)/(L**2)
        self.ratio = v0/R
        self.Fs_x, self.Fs_y, self.Fl_x, self.Fl_y, self.Fs_pwx, self.Fs_pwy = np.zeros(N), np.zeros(N), np.zeros(N), np.zeros(N), np.zeros(N), np.zeros(N)

        orientations = []
        xs = []
        ys = []
        colors = list(matplotlib.colors.CSS4_COLORS.keys())
        colors = [color for color in colors if not self.is_light_color(matplotlib.colors.CSS4_COLORS[color])]
        np.random.shuffle(colors)
        self.colors = colors[:N]

        while len(orientations) < N:
            orientation = np.random.uniform(0, 2*np.pi)
            orientations.append(orientation)
        while len(xs) < N:
            x = np.random.uniform(-(L-2*radius)/2, (L-2*radius)/2)
            y = np.random.uniform(-(L-2*radius)/2, (L-2*radius)/2)
            #Each particle must have a unique initial position
            if not any(np.isclose(x, xs, atol=1.5*radius)) and not any(np.isclose(y, ys, atol=1.5*radius)):
                xs.append(x)
                ys.append(y)

        self.particles = []
        for i in range(len(xs)):
            self.particles.append(Squirmer(xs[i], ys[i], orientations[i], radius, beta, v0))

    def dist_particles(self, particle):
        #Compute the distance between the particle in argument and the other ones
        dist = [self.distance(particle, p) for p in self.particles]
        return dist

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
    
    def forcesSteric(self, particle1, particle2):
        a = self.radius
        Dx = self.vector_x(particle1, particle2)
        Dy = self.vector_y(particle1, particle2)
        dist = self.distance(particle1, particle2)

        tmp = -3*(self.Es/a)*(2*(2*a/dist)**13-(2*a/dist)**7)/np.sqrt(dist)
        print(f"tmp = {tmp}")
        Fs_x =  Dx
        Fs_y = Dy
        return Fs_x, Fs_y

    def ref_border_x(self, particle, boundary):
        #reflective border x
        diff = abs(self.L/2 - abs(particle.x))
        particle.orientation = np.pi - particle.orientation
        if boundary == 1:
            #1 for the right border
            particle.x = self.L/2 - diff
        else:
            particle.x = -self.L/2 + diff

        return particle.x, particle.orientation
    
    def ref_border_y(self, particle, boundary):
        #reflective border y
        particle.orientation = -particle.orientation
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
        for particle in self.particles:
            avrg_orient = self.average_orient(particle)
            noise = np.random.uniform(-self.noise/2, self.noise/2)
            new_orient = np.arctan2(np.sin(avrg_orient), np.cos(avrg_orient)) + noise
            particle.orientation = new_orient

    def compute_forces(self):
        #Compute the steric forces of every particle
        for i, particle in enumerate(self.particles):
            dist = self.dist_particles(particle)
            for j in range(len(dist)):
                if dist[j] != 0 and dist[j]<=self.R:
                    Fs_x, Fs_y = self.forcesSteric(particle, self.particles[j])
                    particle.x -= self.dt*Fs_x
                    particle.y -= self.dt*Fs_y
                    self.particles[j].x += self.dt*Fs_x
                    self.particles[j].y += self.dt*Fs_y        

    def update_position(self):
        #Update position of each particle
        self.compute_forces()
        for i, particle in enumerate(self.particles):
            # if self.Fs_x[i] != 0:
            #     print(f"Fs_x{i} = {self.Fs_x[i]}")
            # if self.Fs_y[i] != 0:
            #     print(f"Fs_y{i} = {self.Fs_y[i]}")
            particle.x += self.v0*self.dt*(np.cos(particle.orientation) + self.Fs_x[i])
            particle.y += self.v0*self.dt*(np.sin(particle.orientation) + self.Fs_y[i])

            if self.L/2 < particle.x + self.radius:
                particle.x, particle.orientation = self.ref_border_x(particle, 1)
            if -self.L/2 > particle.x - self.radius:
                particle.x, particle.orientation = self.ref_border_x(particle, 2)
            if self.L/2 < particle.y + self.radius:
                particle.y, particle.orientation = self.ref_border_y(particle, 1)
            if -self.L/2 > particle.y - self.radius:
                particle.y, particle.orientation = self.ref_border_y(particle, 2)

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
noise = 0.1
Es = 1

vicsek_model = Vicsek_continous(N, R, L, v0, beta, radius, Es, T, dt, noise)

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

#Runs the simulation and plot at intervals
num_steps = 3
for step in range(num_steps):
    start_time = time.time()
    vicsek_model.loop_time()
    end_time = time.time()
    sim_time = end_time - start_time
    print(f"Simulation {step + 1} took {sim_time:.2f} seconds")
    # w = 0
    # for particle in vicsek_model.particles:
    #     print(f"x{w} = {particle.x}")
    #     print(f"y{w} = {particle.y}")
    #     w+=1

    polar = vicsek_model.polar_order_parameter()
    print(f"polar parameter = {polar}")

    fig, ax = plt.subplots(figsize=(8, 8))
    vicsek_model.plot(ax)
    plt.title(f"Positions at Step {step + 1}")
    plt.savefig(f"vicsek_positions_step_{step + 1}.png")
    plt.close()