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
        self.eps = 1e-3
    """
    Desc: Calcule la contribution de l'aimant a un point rpos

    Parametres:
    rpos : Position a laquel on veut la contribution
    RETURN : retourne le champs au point rpos sous la forme d'un numpy array (0,0,Bz)
    """
    def compute_contribution(self, rpos ):

        if distance.euclidean(self.pos,rpos) < self.eps:
            Bz = 0;
        else:
            Bz =  ( Elephant.MU0 / (4*Elephant.PI) ) * ( -self.m ) / ( np.power(distance.euclidean(self.pos, rpos),3) )
            
        return np.array([0., 0., Bz])
