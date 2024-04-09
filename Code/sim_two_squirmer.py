import numpy as np
import matplotlib
import copy

matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

v0 = 1
L = 2
R = L/2
x1, y1 = -0.2, 0.1
x2, y2 = 0.2, 0
orient1, orient2 = 2*np.pi, np.pi
a = 0.05
beta = -7.5
dt = 1e-4
lnEps_cr = 5
Es = 1
T = 50
dt_out = 0.1
Eo = (3/10)*v0/a
ds = 2**(7/6)*a
D = 0

class Squirmer:
    def __init__(self,x,y,orientation,radius,beta,velocity):
        self.x = x
        self.y = y
        self.orientation = orientation
        self.velocity = velocity
        self.radius = radius
        self.beta = beta

    def distance_center(self):
        return np.sqrt(self.x**2 + self.y**2)

    def is_in_square(self,R):
        #return True if the squirmer is in the square
        return (abs(self.x) < (R-self.radius)) and (abs(self.y) < (R-self.radius))

def init_two_squirmers(R,x1,y1,x2,y2,orient1,orient2,radius=a,beta=beta,velocity=v0):

    squirmer1 = Squirmer(x1,y1,orient1,radius,beta,velocity)
    squirmer2 = Squirmer(x2,y2,orient2,radius,beta,velocity)

    if (squirmer1.is_in_square(R) == False) or (squirmer2.is_in_square(R) == False):
        raise ValueError("Squirmers must be inside the square")
    
    square = [R,R]

    return square,squirmer1,squirmer2

def distance_sq(squirmer1, squirmer2):
    Dx = squirmer2.x - squirmer1.x
    Dy = squirmer2.y - squirmer1.y
    return Dx, Dy, np.sqrt(Dx**2 + Dy**2)

def plot_squirmers(R, history):
    plt.figure(figsize=(8, 8))
    plt.plot([-R, R], [-R, -R], 'k-', linewidth=2)  # Bottom side
    plt.plot([-R, R], [R, R], 'k-', linewidth=2)  # Top side
    plt.plot([-R, -R], [-R, R], 'k-', linewidth=2)  # Left side
    plt.plot([R, R], [-R, R], 'k-', linewidth=2)  # Right side

    for step in history:
        plt.scatter(step['squirmer1'].x, step['squirmer1'].y, color='blue', s=100)
        plt.quiver(step['squirmer1'].x, step['squirmer1'].y, np.cos(step['squirmer1'].orientation), np.sin(step['squirmer1'].orientation), color='blue', scale=10)

        plt.scatter(step['squirmer2'].x, step['squirmer2'].y, color='red', s=100)
        plt.quiver(step['squirmer2'].x, step['squirmer2'].y, np.cos(step['squirmer2'].orientation), np.sin(step['squirmer2'].orientation), color='red', scale=10)
    
    plt.axis('equal')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Positions and Orientations of Squirmers')
    plt.grid(True)
    plt.show()


def compute_force_squirmer_border(squirmer, Rpart, R, a, Es):
    tmp = -24*(Es/a)*(2*(a/(R-Rpart))**13-(a/(R-Rpart))**7)/Rpart
    return tmp*squirmer.x, tmp*squirmer.y

def compute_torque_squirmer_border(squirmer, Rpart, R, a, Eo, lnEps_cr):
    ex = squirmer.x / Rpart
    ey = squirmer.y / Rpart

    lnEps = min(lnEps_cr, -np.log((R - Rpart)/a - 1))

    gamma_w = 2 * Eo * (1 + squirmer.beta * (np.cos(squirmer.orientation) * ex + np.sin(squirmer.orientation) * ey)) * \
            lnEps * (np.sin(squirmer.orientation) * ex - np.cos(squirmer.orientation) * ey)
    
    return gamma_w

#Reflective boundary condition
def ref_bound(squirmer, R, a, R2):
    phi_part = np.arctan2(squirmer.y, squirmer.x)
    R_part = np.sqrt(R2)
    x = (2*(R - a) - R_part) * np.cos(phi_part)
    y = (2*(R - a) - R_part) * np.sin(phi_part)
    return x, y


def loop_time(squirmer1, squirmer2, dt, T, R=R, dt_out=dt_out, Es=Es, a=a, Eo=Eo, lnEps_cr=lnEps_cr):
    tout = dt_out
    history = []
    for t in np.arange(dt, T+1, dt):
        Fs_x = 0
        Fs_y = 0
        Dx, Dy, dist = distance_sq(squirmer1, squirmer2)
        #print("dist =", dist)
        #print("ds =", ds)
        #Force between squirmers
        if dist < ds:
            tmp = -12*(Es/a)*(2*(2*a/dist)**13-(2*a/dist)**7)/dist
            #print("tmp =", tmp)
            Fs_x = tmp * Dx
            #print("Fs_x =", tmp)
            Fs_y = tmp * Dy

        #Force between a squirmer and a border
        Rpart1 = squirmer1.distance_center()
        Rpart2 = squirmer2.distance_center()
        Fs_pw1 = [0,0]
        Fs_pw2 = [0,0]
        if ((R-Rpart1) < 2**(1/6)*a):
            Fs_pw1[0], Fs_pw1[1] = compute_force_squirmer_border(squirmer1, Rpart1, R, a, Es)
        if ((R-Rpart2) < 2**(1/6)*a):
            Fs_pw2[0], Fs_pw2[1] = compute_force_squirmer_border(squirmer2, Rpart2, R, a, Es)
        
        #Compute torques exerted on squirmer by other squirmer
        val1 = 0
        val2 = 0
        if dist < 3*a:
            ex = Dx/dist
            ey = Dy/dist

            lnEps = min(lnEps_cr, -np.log(dist/a - 2))
            val1 = Eo * (1 + beta * (np.cos(squirmer1.orientation) * ex + np.sin(squirmer1.orientation) * ey)) * \
                    lnEps * (ex * np.sin(squirmer1.orientation) - ey * np.cos(squirmer1.orientation))
            val2 = Eo * (1 + beta * (np.cos(squirmer2.orientation) * ex + np.sin(squirmer2.orientation) * ey)) * \
                    lnEps * (ex * np.sin(squirmer2.orientation) - ey * np.cos(squirmer2.orientation))

        #Compute torque exerted on squirmer by the wall
        gamma_w1 = 0
        gamma_w2 = 0
        if ((R - Rpart1) < 2**(1/6) * a):
            gamma_w1 = compute_torque_squirmer_border(squirmer1, Rpart1, R, a, Eo, lnEps_cr)
        if ((R - Rpart2) < 2**(1/6) * a):
            gamma_w2 = compute_torque_squirmer_border(squirmer2, Rpart2, R, a, Eo, lnEps_cr)
    
        #Evolution of position
        squirmer1.x += dt*(squirmer1.velocity * np.cos(squirmer1.orientation) + Fs_x + Fs_pw1[0])
        squirmer1.y += dt*(squirmer1.velocity * np.sin(squirmer1.orientation) + Fs_y + Fs_pw1[1])
        squirmer1.orientation += dt*(val1 + 0.25*val1 + gamma_w1)

        squirmer2.x += dt*(squirmer2.velocity * np.cos(squirmer2.orientation) - Fs_x + Fs_pw2[0])
        squirmer2.y += dt*(squirmer2.velocity * np.sin(squirmer2.orientation) - Fs_y + Fs_pw2[1])
        squirmer2.orientation += dt*(val2 + 0.25*val2 + gamma_w2)

        #Reflective boundary
        R2_sq1 = squirmer1.x**2 + squirmer1.y**2
        #print("x et y:",squirmer1.x, squirmer1.y,"\n")
        #print("R2", R2_sq1,"\n")
        if R2_sq1 > (R-a)**2:
            squirmer1.x, squirmer1.y = ref_bound(squirmer1, R, a, R2_sq1)
            #print(squirmer1.x, squirmer1.y,"\n")
        R2_sq2 = squirmer2.x**2 + squirmer2.y**2
        if R2_sq2 > (R-a)**2:
            squirmer2.x, squirmer2.y = ref_bound(squirmer2, R, a, R2_sq2)
            #print(squirmer2.x, squirmer2.y,"\n")

        #Plots
        if t >= tout:
            sq1_copie = copy.deepcopy(squirmer1)
            sq2_copie = copy.deepcopy(squirmer2)
            history.append({'squirmer1':sq1_copie, 'squirmer2':sq2_copie})
            tout += dt_out
    return history

square, squirmer1, squirmer2 = init_two_squirmers(R, x1, y1, x2, y2, orient1, orient2)
history = loop_time(squirmer1, squirmer2, dt, T, R)
plot_squirmers(R, history)