import numpy as np
import matplotlib.pyplot as plt

class Squirmer:
    def __init__(self,x,y,orientation,radius,beta):
        self.x = x
        self.y = y
        self.orientation = orientation
        self.radius = radius
        self.beta = beta

    def is_in_square(self,L):
        #return True if the squirmer is in the square
        return (abs(self.x) < (L/2-self.radius)) and (abs(self.y) < (L/2-self.radius))

def init_two_squirmers(L,x1,y1,x2,y2,orient1,orient2,radius=0.05,beta1=-7.5,beta2=-7.5):

    squirmer1 = Squirmer(x1,y1,orient1,radius,beta1)
    squirmer2 = Squirmer(x2,y2,orient2,radius,beta2)

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

def loop_time(square, squirmer1, squirmer2):
    dist = distance_sq(squirmer1, squirmer2)

    return 0

L = 2
x1, y1 = -0.5, 0.5
x2, y2 = 0.5, -0.5
orient1, orient2 = np.pi/4, np.pi/3

square, squirmer1, squirmer2 = init_two_squirmers(L, x1, y1, x2, y2, orient1, orient2)
plot_squirmers(square, squirmer1, squirmer2)