from Code.squirmer import Squirmer
import numpy as np

def test_sq():
    x, y = 0.5, 0.2
    orient = np.pi
    a = 0.05
    beta = -2
    v0 = 0.7
    squirmer = Squirmer(x,y,orient,a,beta,v0)
    assert squirmer.B1 == (3./2.)*v0
    assert squirmer.B2 == beta*(3./2.)*v0