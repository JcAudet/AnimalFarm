import numpy as np
from Panda import Panda
from Poulpe2 import Poulpe

from numpy.linalg import norm

## classe objet du solveur

class Dauphin:

    """
    Desc: Solveur Verlet

    Parametres:
        Panda : Puck bougeants
        Poulpe : Outils pour le magnetic field
        dt : discretisation en temps
        dq : discretisation de l'espace
        xlim, ylim : np.array([min,max])
    Returns :
        positions, velocities : np.array([N][x,y])
    """
    def __init__(self, panda, poulpe, dt, dq, xlim, ylim):
        self.panda = panda
        self.poulpe = poulpe
        self.dt = dt # time step
        self.dq = dq # space step
        self.xmin, self.xmax = xlim[0], xlim[1]
        self.ymin, self.ymax = ylim[0], ylim[1]

    def solve(self, N):
        """
        Solve for N time steps
        """
        positions = np.zeros([N,2])
        velocities = np.zeros([N,2])
        n=0
        positions[n,:] = self.panda.pos0
        velocities[n,:] = self.panda.vit0
        # tn=1 step -------------------------------------------------
        n=1
        pos1, vit1 = verlet_step_1(self.panda, self.poulpe, self.dt, self.dq, self.xmin, self.xmax, self.ymin, self.ymax)
        positions[n,:] = pos1
        velocities[n,:] = vit1
        # n>1 steps ------------------------------------------------
        if N > 1:
            for n in range(N-1):
                n +=1
                posNplus1, vitNplus1 = verlet_step_n(self.panda, self.poulpe, self.dt, self.dq, self.xmin, self.xmax, self.ymin, self.ymax)
                positions[n,:] = posNplus1
                velocities[n,:] = vitNplus1

        return positions, velocities


def puck_outside(posNplus1, xmin, xmax, ymin, ymax):
    """ Check if the position of the puck is outside the limits.
    If is outside, return True """
    if  posNplus1[0]<xmin or posNplus1[0]>xmax or ymin>posNplus1[1] or posNplus1[0]>ymax:
        return True

def calculate_force(panda, pos, poulpe, dq):
    """Calculates the energy derivatives in x and y in space at the position pos. FORCE!!
    Derivative :
            U(q) = -m@B(q) -> U(q+dq)=-m@B(q+dq)
            dU/dq_i = (-1)*( m@B(q+dq_i) - m@B(q - dq_i) )/(2*dq_i)
    Args:
        pos : np.array([qx,qy])
        dq : float
        m : np.array([mx, my, mz])
    Returns:
        vecDerivative : np.array([dB/dx, dB/dy])
                        Derivatives of B in direction x and y at point pos.
    """
    vecdqX = np.array([dq,0])
    vecdqY = np.array([0,dq])
    derivativeX = ( poulpe.compute_field(pos+vecdqX)@panda.m -poulpe.compute_field(pos-vecdqX)@panda.m )/(2*dq)
    derivativeY = ( poulpe.compute_field(pos+vecdqY)@panda.m -poulpe.compute_field(pos-vecdqY)@panda.m)/(2*dq)
    vecDerivative = (-1)*np.array([derivativeX, derivativeY])
    return vecDerivative


def verlet_step_1(panda, poulpe, dt, dq, xmin, xmax, ymin, ymax):
    """First step of Verlet integration
    Calculates position at time n=1 from intial conditions

    q1 = q0 + v0*dt + 0.5*dt^2*dU(q0)/dq

    Args:
        panda (Panda)
        poulpe (Poulpe)
        dt (float) : time step
        dq (float) : space step
    Returns:
        pos1 : np.array([qx, qy])
    """
    force0=calculate_force(panda,panda.pos, poulpe, dq)
    pos1 = panda.pos + panda.vit*dt +(0.5)*(dt**2)*force0/panda.mass
    # Update position
    panda.update_pos(pos1)
    # Calcul vitesse
    force1 = calculate_force(panda, panda.pos, poulpe, dq)
    vit1 = panda.vit + 0.5*dt*(force0 + force1)/panda.mass
    # Update velocity
    panda.update_vit(vit1)
    # outside ?
    # if puck_outside(pos1, xmin, xmax, ymin, ymax):
    #     pass
    return panda.pos, panda.vit

def verlet_step_n(panda, poulpe, dt, dq, xmin, xmax, ymin, ymax):
    """ Step n>1 of Verlet integration
    Calculates next position from current position

    qn+1 = 2*qn - qn-1 + dt^2*dU(qn)/dq

    Args:
        panda (Panda)
        poulpe (Poulpe)
        dt (float) : time step
        dq (float) : space step
    Returns:
        posn : np.array([qx, qy])
    """
    forceN=calculate_force(panda, panda.pos, poulpe, dq)
    posNplus1 = 2*panda.pos - panda.lastPos + (dt**2)*forceN/panda.mass
    # Update position
    panda.update_pos(posNplus1)
    # calculer vitesse
    forceNplus1 = calculate_force(panda, panda.pos, poulpe, dq)
    vitNplus1 = panda.vit + 0.5*dt*(forceN + forceNplus1)/panda.mass
    # Update velocity
    panda.update_vit(vitNplus1)
    # outside ?
    # if puck_outside(posNplus1, xmin, xmax, ymin, ymax):
    #     pass
    return panda.pos, panda.vit
