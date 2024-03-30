import numpy as np

class Squirmer:
    def __init__(self,x,y,orientation,radius,beta):
        self.x = x
        self.y = y
        self.orientation = orientation
        self.radius = radius
        self.beta = beta

    def is_in_square(self,L):
        #return True if the squirmer is in the square
        return (abs(self.x) <= (L-self.radius)) and (abs(self.y) <= (L-self.radius))

def init_two_squirmers(L,x1,y1,x2,y2,orient1,orient2,radius=0.05,beta1=-7.5,beta2=-7.5):

    squirmer1 = Squirmer(x1,y1,orient1,radius,beta1)
    squirmer2 = Squirmer(x2,y2,orient2,radius,beta2)

    if (squirmer1.is_in_square(L) == False) or (squirmer2.is_in_square(L) == False):
        raise ValueError("Squirmers must be inside the square")
    
    square = [L,L]

    return square,squirmer1,squirmer2
