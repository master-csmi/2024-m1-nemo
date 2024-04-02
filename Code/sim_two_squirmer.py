import numpy as np
import matplotlib
# Utiliser un backend compatible avec l'affichage graphique
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

v0 = 1
L = 2
x1, y1 = -0.5, 0.5
x2, y2 = 0.5, -0.5
orient1, orient2 = np.pi/4, np.pi/3
a = 0.05
beta = -7.5
dt = 1e-4
lnEps_cr = 5
Es = 1
T = 2
dt_out = 0.1
Eo = (3/10)*v0/a
ds = 2^(7/6)*a

class Squirmer:
    def __init__(self,x,y,orientation,radius,beta,velocity=v0):
        self.x = x
        self.y = y
        self.orientation = orientation
        self.velocity = velocity
        self.radius = radius
        self.beta = beta

    def distance_center(self):
        return np.sqrt(self.x^2 + self.y^2)

    def is_in_square(self,L):
        #return True if the squirmer is in the square
        return (abs(self.x) < (L/2-self.radius)) and (abs(self.y) < (L/2-self.radius))

def init_two_squirmers(L,x1,y1,x2,y2,orient1,orient2,radius=a,beta1=beta,beta2=beta,velocity=v0):

    squirmer1 = Squirmer(x1,y1,orient1,radius,beta1,velocity)
    squirmer2 = Squirmer(x2,y2,orient2,radius,beta2,velocity)

    if (squirmer1.is_in_square(L) == False) or (squirmer2.is_in_square(L) == False):
        raise ValueError("Squirmers must be inside the square")
    
    square = [L,L]

    return square,squirmer1,squirmer2

def distance_sq(squirmer1, squirmer2):
    Dx = squirmer2.x - squirmer1.x
    Dy = squirmer2.y - squirmer2.y
    return np.square(Dx**2 + Dy**2)

#def steric_force(squirmer1, squirmer2):
    #distance = distance_sq(squirmer1, squirmer2)
    #if distance

def plot_squirmers(square, squirmer1, squirmer2):
    plt.figure(figsize=(8, 8))
    plt.plot([-square[0]/2, square[0]/2], [-square[1]/2, -square[1]/2], 'k-', linewidth=2)  # Bottom side
    plt.plot([-square[0]/2, square[0]/2], [square[1]/2, square[1]/2], 'k-', linewidth=2)  # Top side
    plt.plot([-square[0]/2, -square[0]/2], [-square[1]/2, square[1]/2], 'k-', linewidth=2)  # Left side
    plt.plot([square[0]/2, square[0]/2], [-square[1]/2, square[1]/2], 'k-', linewidth=2)  # Right side

    # Plot squirmer 1
    plt.scatter(squirmer1.x, squirmer1.y, color='blue', s=100)
    plt.quiver(squirmer1.x, squirmer1.y, np.cos(squirmer1.orientation), np.sin(squirmer1.orientation), color='blue', scale=10)

    # Plot squirmer 2
    plt.scatter(squirmer2.x, squirmer2.y, color='red', s=100)
    plt.quiver(squirmer2.x, squirmer2.y, np.cos(squirmer2.orientation), np.sin(squirmer2.orientation), color='red', scale=10)

    plt.axis('equal')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Positions and Orientations of Squirmers')
    plt.grid(True)
    plt.show()

def loop_time(square, squirmer1, squirmer2, dt, T):
    for t in np.arange(0, T, dt):
        Fs_x = 0
        Fs_y = 0
        dist = distance_sq(squirmer1, squirmer2)
        #Force between squirmers
        if dist < ds:
            tmp = -12*(Es/a)*(2*(2*a/dist)^13-(2*a/dist)^7)/dist
            Fs_x = tmp * (squirmer2.x - squirmer1.x)
            Fs_y = tmp * (squirmer2.y - squirmer1.y)

        #Force between a squirmer and a border
        Rpart1 = squirmer1.distance_center()
        Rpart2 = squirmer2.distance_center()
        Fs_pw1 = []
        Fs_pw2 = []
        if ((L-Rpart1) < 2^(1/6)*a):
            tmp = -24*(Es/a)*(2*(a/(L-Rpart1))^13-(a/(L-Rpart1))^7)/Rpart1
            Fs_pw1.append(tmp*squirmer1.x)
            Fs_pw1.append(tmp*squirmer1.y)
        if ((L-Rpart2) < 2^(1/6)*a):
            tmp = -24*(Es/a)*(2*(a/(L-Rpart2))^13-(a/(L-Rpart2))^7)/Rpart2
            Fs_pw2.append(tmp*squirmer1.x)
            Fs_pw2.append(tmp*squirmer1.y)
    return 0

square, squirmer1, squirmer2 = init_two_squirmers(L, x1, y1, x2, y2, orient1, orient2)
plot_squirmers(square, squirmer1, squirmer2)