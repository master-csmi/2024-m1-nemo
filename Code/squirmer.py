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