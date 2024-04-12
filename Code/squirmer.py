import numpy as np

class Squirmer:
    def __init__(self,x,y,orientation,radius,beta,velocity):
        self.x = x
        self.y = y
        self.orientation = orientation
        self.velocity = velocity
        self.radius = radius
        self.beta = beta
        self.B1 = (2/3)*velocity
        self.B2 = beta*(2/3)*velocity

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