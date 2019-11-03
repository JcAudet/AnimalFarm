import numpy as np
import matplotlib.tri as tri
from Calmar2 import Calmar

## Classe utiliser pour obtenir le champs magnetique en tout point grace a une liste de aimants (Calmar)

############## SSSSSSSSQQQQQQQQQQQQUUUUUUUUUIIIIIIIIIIIIDDDDDDDDDDDDD ###################

class Poulpe:

    """
    Desc: Constructeur
    Parametres:
        vision : object graphique du GUI pour la fonction update_Mesh
    """
    def __init__(self, vision=None):
        self.listeCalmar = []
        self.vision = vision

    """
    Desc : Ajoute un aimant a la liste d'aimant

    Parametres:
        pos : position de l'aimant a ajouter
        m : Aimantation de l'aimant a ajouter a la table
    """
    def add_calmar(self, pos, m):

        self.listeCalmar.append( Calmar(pos, m) )

    """
    Desc : Ajoute plusieurs calmar d'un seul cout a la liste de Calmar

    Parametres:
        calmars : liste de Calmars a ajouter a SSQQUUUIIIDDD
    """
    def add_calmar_liste(self, calmars):

        for calmar in calmars:
            self.listeCalmar.append(calmar)

    """
    Desc : remplace la liste de Calmars presente par une nouvelle liste de Calmars

    Parametres:
        calmars : nouvelle liste de calmars
    """
    def set_calmars(self, calmars):

        self.listeCalmar = calmars

        self.update_mesh()

    """
    Desc : Calcule le champs totale de tous les aimant a un point pos

    Parametres:
        pos : position du champs desire
        RETURN : retourne un numpy array du champs total de tous les aimants
    """
    def compute_field(self, pos ):

        Btot = np.array([0.,0.,0.])

        for calmar in self.listeCalmar:

            Btot = Btot + calmar.compute_contribution(pos)

        
        return Btot

    """
    Desc : Update values dans le graph GUI VISION grace a la fonction compute_field
    """
    def update_mesh(self):

        x = self.vision.axes.get_xlim()
        y = self.vision.axes.get_ylim()
        X = np.linspace( x[0], x[1], 10000 )
        Y = np.linspace( y[0], y[1], 10000 )

        mesh = tri.Triangulation( X, Y )

        i=0
        field = np.zeros(X.len() * Y.len())
        for x in X:
            for y in Y:
                field[i] = self.compute_field(np.array([x,y]))
                i += 1

        cont = self.vision.axes.tricontourf(mesh, field)

        self.vision.update_graph()
