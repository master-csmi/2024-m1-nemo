import numpy as np
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

        orientations = set()
        xs = set()
        ys = set()

        while len(orientations) < N:
            orientation = np.random.uniform(0, 2*np.pi)
            orientations.add(orientation)
        while len(xs) < N:
            x = np.random.uniform(-(L-2*radius)/2, (L-2*radius)/2)
            xs.add(x)
        while len(ys) < N:
            y = np.random.uniform(-(L-2*radius)/2, (L-2*radius)/2)
            ys.add(y)

        self.particles = set()
        for i in range(len(xs)):
            self.particles.add(Squirmer(xs[i], ys[i], orientations[i], radius, beta, v0))

    def ref_border_x(self, particle, boundary):
        #reflective border x
        diff = self.L/2 - abs(particle.x)
        particle.orientation = np.pi - particle.orientation
        if boundary == 1:
            #1 for the right border
            particle.x = self.L/2 - diff
        else:
            particle.x = self.L/2 + diff
    
    def ref_border_y(self, particle, boundary):
        #reflective border y
        particle.orientation = -particle.orientation
        diff = self.L/2 - abs(particle.y)
        if boundary == 1:
            #1 for the up boundary
            particle.y = self.L/2 - diff
        else:
            particle.y = -self.L/2 + diff

    def distance(self, p1, p2):
        #Compute the distance between two particles
        return np.sqrt((p1.x-p2.x)**2 + (p1.y-p2.y)**2)
    
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
            particle.x += self.v0 * self.dt * np.cos(particle.orientation)
            particle.y += self.v0 * self.dt * np.sin(particle.orientation)

            if particle.x >= self.L/2:
                self.ref_border_x(particle, 1)
            if particle.x <= -self.L/2:
                self.ref_border_x(particle, 2)
            if particle.y >= self.L/2:
                self.ref_border_x(particle, 1)
            if particle.y <= -self.L/2:
                self.ref_border_x(particle, 2)

    def loop_time(self):
        for _ in np.arange(self.dt, self.T, self.dt):
            self.update_orient()
            self.update_position()
