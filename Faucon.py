from Poulpe import *
from Panda import *
from Dauphin import *
from OurangOutan import *

## Classe agissant comme Overhead de l'algorithme de resolution

class Faucon:

	#Valeurs par défault (initiales) des paramètres numériques


	# Constructeur
	def __init__(self, lionceau, dictHyene, course, vision, fourmi):
		self.nt = fourmi.nt
		self.dq = fourmi.dq
		self.dt = fourmi.dt
		self.poulpe = Poulpe(vision)
		self.dictHyene = dictHyene
		self.dict_hyene_to_poulpe()
		#self.panda = Panda(lionceau.pos, lionceau.vit, lionceau.force,lionceau.masse)
		self.dauphin = Dauphin(self.panda, self.poulpe, fourmi.dt,fourmi.dq, vision.axes.get_xlim(), vision.axes.get_ylim())

	def run_simulation(self, pos0, vit0, m0):
                
                
                panda = Panda(pos0, vit0, m0)

		xlim = np.array([-100,100])
		ylim = np.array([-100,100])

		dauphin = Dauphin(panda, self.poulpe, dt, dq, xlim, ylim)

		sol = dauphin.solve(self.nt)

		return sol


#Transforme un dictionnaire de Hyene en liste de Calmar et le set comme liste de calmars de poulpe
	def dict_hyene_to_poulpe(self):
		calmarList = []
		for key in self.dictHyene:
			calmarList.append(Calmar(self.dictHyene[key].force,self.dictHyene[key].pos))
		self.poulpe.set_calmars(calmarList)

	#Mise a jour du graphique
	def update_graphB(self):
		self.poulpe.update_mesh()
