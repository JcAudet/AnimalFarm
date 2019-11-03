import numpy as np
import matplotlib.tri as tri
from Calmar import Calmar
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.colors import SymLogNorm
from matplotlib import ticker, cm
plt.ion()
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
        self.vision.axes = self.vision.Fig.gca()
        self.vision.update_graph()
        self.initiate_mesh()

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

            Btot += calmar.compute_contribution(pos)

        return Btot

    def initiate_mesh(self):


        x = self.vision.axes.get_xlim()
        y = self.vision.axes.get_ylim()
        self.X = np.linspace( x[0], x[1], 150 )
        self.Y = np.linspace( y[0], y[1], 150 )

        self.xx, self.yy = np.meshgrid(self.X,self.Y)

        #mesh = tri.Triangulation( X, Y )

        #print(mesh)

        i=0
        field = np.zeros((self.X.size, self.Y.size))
        for x in self.X:
            j=0
            for y in self.Y:
                field[i][j] = self.compute_field(np.array([x,y]))[2]
                j += 1
            i+=1


        print('init')
        #cont = self.vision.axes.tricontourf(mesh, field)

        self.vision.axes.pcolormesh(self.xx, self.yy,field,cmap='RdBu',
                                    norm = colors.SymLogNorm(linthresh=1e-6,
                                                             linscale = 10, vmin=-np.abs(field).max(),vmax=np.abs(field).max()) )


        #self.img = self.vision.axes.imshow(field, cmap='RdBu',
        #                                   aspect='auto',
        #                                   interpolation='none',
        #                                   norm = SymLogNorm(linthresh=1e-8,
        #                                                     linscale = 10,
        #                                                     vmin=-np.abs(field).max(),
        #                                                     vmax=np.abs(field).max()),
        #                                  extent=(self.X[0], self.X[-1],
        #                                          self.Y[0], self.Y[-1]),
        #                                   )
#        self.Quad = self.vision.axes.pcolormesh(self.xx,
#                                                self.yy,field,cmap='RdBu',
#                                                shading = 'gouraud',
#                                                norm = colors.SymLogNorm(linthresh=1e-6,
#                                                                         linscale = 5,
#                                                                         vmin=-np.abs(field).max(),
#                                                                         vmax=np.abs(field).max()),
#                                               animated=True,
#                                               figure=self.vision.Fig
#                                               )
#
        self.vision.update_graph()
        #self.Quad.draw()

    """
    Desc : Update values dans le graph GUI VISION grace a la fonction compute_field
    """
    def update_mesh(self):


        i=0
        field = np.zeros((self.X.size, self.Y.size))
        for x in self.X:
            j=0
            for y in self.Y:
                field[i][j] = self.compute_field(np.array([x,y]))[2]
                j += 1
            i+=1
        #self.img.set_data(field)
        print('update')
        #self.Quad.set_array(field.ravel())
        #cont = self.vision.axes.tricontourf(mesh, field)
        self.vision.axes.pcolormesh(self.xx, self.yy,field,cmap='RdBu',
                                    shading = 'gouraud',
                                    norm = colors.SymLogNorm(linthresh=1e-6,
                                                            #linscale = 10, vmin=-np.abs(field).max(),vmax=np.abs(field).max()) )
                                                             linscale = 5,
                                                             vmin=-np.abs(field).max(),
                                                             vmax=np.abs(field).max()),
                                   )

        self.vision.update_graph()
