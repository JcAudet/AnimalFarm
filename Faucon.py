from Poulpe import *
from Panda import *
from Dauphin import *
from OurangOutan import *

## Classe agissant comme Overhead de l'algorithme de resolution

class Faucon:

	#Valeurs par défault (initiales) des paramètres numériques


	# Constructeur
	def __init__(self, lionceau,dictHyene,course, vision,fourmi):
		self.nt = fourmi.nt
		self.poulpe = Poulpe(vision)
		self.dictHyene = dictHyene
		self.dict_hyene_to_poulpe()
		self.panda = Panda(lionceau.pos, lionceau.vit, lionceau.force,lionceau.masse)
		self.dauphin = Dauphin(self.panda, self.poulpe, fourmi.dt,fourmi.dq, vision.axes.get_xlim(), vision.axes.get_ylim())

	def run_simulation(self):
		self.dauphin.solve(self.nt.get())


#Transforme un dictionnaire de Hyene en liste de Calmar et le set comme liste de calmars de poulpe
	def dict_hyene_to_poulpe(self):
		calmarList = []
		for key in self.dictHyene:
			calmarList.append(Calmar(self.dictHyene[key].force,self.dictHyene[key].pos))
		self.poulpe.set_calmars(calmarList)

	#Mise a jour du graphique
	def update_graphB(self):
		self.poulpe.update_mesh()
