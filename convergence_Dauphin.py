from Calmar2 import Calmar
import Elephant
from Dauphin import Dauphin
from Panda import Panda
from Poulpe2 import Poulpe

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit


# POOOOUUUUUUUUUULLLLLPPPPPPPPPPEEEEEEEEEEEE
poulpe = Poulpe()
poulpe.add_calmar(10, np.array([0,1e-1]))
# poulpe.add_calmar(-10, np.array([0.5,0.5]))
# poulpe.add_calmar(50, np.array([1,1]))

# PPPAAAAAAANNNNNNNNDDDDDDDAAAAAAAAAAA
pos0 = np.array([0,0])
vit0 = np.array([1e-5,0])
m = np.array([0,0,1])
panda = Panda(pos0, vit0, m)
panda2 = Panda(pos0, vit0, m)

# # Convergence Temporelle
dq = 1e-1
Lt = 1
Nts = np.power(2, np.arange(3,13)) #Differents nombre de points
dts = np.zeros(len(Nts)) #dt pour les differents nombre de points
errsdt = np.zeros(len(Nts)) #erreurs pour les differents nombre de points
errsvit = np.zeros(len(Nts))
xlim = np.array([-100,100])
ylim = np.array([-100,100])

for indexdt in range(len(Nts)): # index sur les Nt
    Nt = Nts[indexdt] #nbre de points de temps
    Nt2 = Nt*2
    dt = Lt/Nt
    dt2 = Lt/Nt2
    dts[indexdt] = dt

    dauphin = Dauphin(panda, poulpe, dt, dq, xlim, ylim)
    dauphin2 = Dauphin(panda2, poulpe, dt2, dq, xlim, ylim)

    dauphin.solve(Nt)
    dauphin2.solve(Nt2)

    sol = np.array(dauphin.panda.storePos) # sol[temps][x,y]
    sol2 = np.array(dauphin2.panda.storePos)

    vit = np.array(dauphin.panda.storeVit) # sol[temps][x,y]
    vit2 = np.array(dauphin2.panda.storeVit)

    dauphin.panda.clear_panda()
    dauphin2.panda.clear_panda()

    errxy=0
    errvit=0
    for t in range(Nt): #indice sur le temps
        errxy += (np.linalg.norm(sol[t] - sol2[2*t]))/Nt
        errvit += (np.linalg.norm(vit[t] - vit2[2*t]))/Nt
    errsdt[indexdt] = errxy
    errsvit[indexdt] = errvit


def f(x,p,c):
    return c + p*np.log(x)

popt, pcov = curve_fit(f, dts, np.log(errsdt))
print(popt)
poptV, pcovV = curve_fit(f, dts, np.log(errsvit))
print(poptV)

plt.figure()
plt.title('Analyse de convergence temporelle')
plt.xlabel('dt')
plt.ylabel('Erreur')
plt.loglog(dts, np.exp(f(dts, popt[0], popt[1])))
plt.loglog(dts, np.exp(f(dts, poptV[0], poptV[1])))
plt.loglog(dts, errsdt, 'o')
plt.loglog(dts, errsvit, 'x')
plt.show()
plt.figure()
plt.title('Trajectoire (Analyse spatiale)')
plt.xlabel('x')
plt.ylabel('y')
plt.plot(sol2[:,0], sol2[:,1], 'o')
plt.show()



# Convergence Spatiale (attention aux PPPAAAAAAANNNNNNNNDDDDDDDAAAAAAAAAAAs)
Lt=1
Nt = int(2**10)
dt = Lt/Nt
dqs = np.array([2**(-11), 2**(-10), 2**(-9), 2**(-8), 2**(-7), 2**(-6),2**(-5)]) #Differents nombre de points

errsdq = np.zeros(len(dqs)) #erreurs pour les differents nombre de points
errsvit = np.zeros(len(dqs))
# print(len(dqs))

for indexdq in range(len(dqs)): # index sur les dq

    dq = dqs[indexdq]
    dq2 = dqs[indexdq]/2

    dauphin = Dauphin(panda, poulpe, dt, dq,xlim,ylim)
    dauphin2 = Dauphin(panda2, poulpe, dt, dq2,xlim,ylim)

    dauphin.solve(Nt)
    dauphin2.solve(Nt)

    sol = np.array(dauphin.panda.storePos) # sol[temps][x,y]
    sol2 = np.array(dauphin2.panda.storePos)

    vit = np.array(dauphin.panda.storeVit) # sol[temps][x,y]
    vit2 = np.array(dauphin2.panda.storeVit)

    dauphin.panda.clear_panda()
    dauphin2.panda.clear_panda()

    errxy=0
    errvit=0
    for i in range(len(dqs)): #indice sur le temps
        errxy += (np.linalg.norm(sol[i] - sol2[i]))/Nt
        errvit += (np.linalg.norm(vit[i] - vit2[i]))/Nt
    errsdq[indexdq] = errxy
    errsvit[indexdq] = errvit


def f(x,p,c):
    return c + p*np.log(x)

popt, pcov = curve_fit(f, dqs, np.log(errsdq))
print(popt)
poptV, pcovV = curve_fit(f, dqs, np.log(errsvit))
print(poptV)

plt.figure()
plt.title('Analyse de convergence spatiale')
plt.xlabel('dq')
plt.ylabel('Erreur')
plt.loglog(dqs, np.exp(f(dqs, popt[0], popt[1])))
plt.loglog(dqs, errsdq, 'o')
plt.loglog(dqs, np.exp(f(dqs, poptV[0], poptV[1])))
plt.loglog(dqs, errsvit, 'x')
plt.show()

plt.figure()
plt.title('Trajectoire (Analyse spatiale)')
plt.xlabel('x')
plt.ylabel('y')
plt.plot(sol2[:,0], sol2[:,1], 'o')
plt.show()
