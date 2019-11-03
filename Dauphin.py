import numpy as np
from Panda import Panda
from Poulpe import Poulpe

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
        # n=1 step -------------------------------------------------
        # Calculate next position and position
        pos1 = verlet_step_1(self.panda, self.poulpe, self.dt, self.dq, self.xmin, self.xmax, self.ymin, self.ymax)
        # n>1 steps ------------------------------------------------
        if N > 1:
            for n in range(N-1):
                # Calculate next position (and velocity - TO DO)
                posn = verlet_step_n(self.panda, self.poulpe, self.dt, self.dq, self.xmin, self.xmax, self.ymin, self.ymax)

def puck_outside(posNplus1, xmin, xmax, ymin, ymax):
    """ Check if the position of the puck is outside the limits.
    If is outside, return True """
    if  posNplus1[0]<xmin or posNplus1[0]>xmax or ymin>posNplus1[1] or posNplus1[0]>ymax:
        return True

def reflection(panda, qoutside, xmin, xmax, ymin, ymax):
    """ Changes posNplus1 and vitNplus1 for reflection """
    pass
    # x = 1 #qnplus1 outside bounds
    # y = 1
    # vx = 1 #vnplus1 outside bounds
    # vy = 1
    #
    # if x>=xmax:
    #     xgood = 2*xmax - x
    #     vxgood = -vx
    # elif x<=xmin:
    #     xgood = 2*xmin - x
    #     vxgood = -vx
    # if y>=ymax:
    #     self.pos[1] = 2*ymax - y
    #     self.vit[1] = -vy
    # elif y<=ymin:
    #     self.pos[1] = 2*ymin - y
    #     self.vit[1] = -vy



def space_derivative_energy(panda, poulpe, dq):
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
    derivativeX = ( poulpe.compute_field(panda.pos+vecdqX)@panda.m -poulpe.compute_field(panda.pos-vecdqX)@panda.m )/(2*dq)
    derivativeY = ( poulpe.compute_field(panda.pos+vecdqY)@panda.m -poulpe.compute_field(panda.pos-vecdqY)@panda.m)/(2*dq)
    vecDerivative = (-1)*np.array([derivativeX, derivativeY])
    return vecDerivative

def forceOutsideBounds(panda, pos, poulpe, dq):
    """Calculates the energy derivatives in x and y in space at the position pos OUTSIDE BOUNDS.
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
    force0=space_derivative_energy(panda, poulpe, dq)
    pos1 = panda.pos + panda.vit*dt +(0.5)*(dt**2)*force0/panda.mass
    if puck_outside(pos1, xmin, xmax, ymin, ymax):
        reflection(panda, pos1, xmin, xmax, ymin, ymax)
    else:
        # Update position
        panda.update_pos(pos1)
        # Calcul vitesse
        force1 = space_derivative_energy(panda, poulpe, dq)
        vit1 = panda.vit + 0.5*dt*(force0 + force1)/panda.mass
        # Update velocity
        panda.update_vit(vit1)
        # outside ?
    return panda.pos

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
    forceN=space_derivative_energy(panda, poulpe, dq)
    posNplus1 = 2*panda.pos - panda.lastPos + (dt**2)*forceN/panda.mass

    # outside ?
    if puck_outside(posNplus1, xmin, xmax, ymin, ymax):
        reflection(panda, posNplus1, xmin, xmax, ymin, ymax)
    else:
        # Update position
        panda.update_pos(posNplus1)
        # calculer vitesse
        forceNplus1 = space_derivative_energy(panda, poulpe, dq)
        vitNplus1 = panda.vit + 0.5*dt*(forceN + forceNplus1)/panda.mass
        # Update velocity
        panda.update_vit(vitNplus1)
    return panda.pos
