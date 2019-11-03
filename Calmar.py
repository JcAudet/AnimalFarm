import numpy as np
import scipy.spatial.distance as distance
import Elephant

## Classe d'objet repressantant les aimants fixe sur la table de jeu

class Calmar:

    """
    Desc: Constructeur

    Parametres:
    pos : Position
    m : force aimantation
    """
    def __init__(self, m, pos ):
        self.pos = pos
        self.m = m
        self.eps = 1e-6
    """
    Desc: Calcule la contribution de l'aimant a un point rpos

    Parametres:
    rpos : Position a laquel on veut la contribution
    RETURN : retourne le champs au point rpos sous la forme d'un numpy array (0,0,Bz)
    """
    def compute_contribution(self, rpos ):

        if distance.euclidean(np.array([self.pos[0].get(), self.pos[1].get()]),rpos) < self.eps:
            Bz =  ( Elephant.MU0 / (4*Elephant.PI) ) * ( -self.m.get() ) / ( np.power(self.eps,3) )
        else:
            Bz =  ( Elephant.MU0 / (4*Elephant.PI) ) * ( -self.m.get() ) / ( np.power(distance.euclidean(np.array([self.pos[0].get(), self.pos[1].get()]), rpos),3) )

        return np.array([0., 0., Bz])
