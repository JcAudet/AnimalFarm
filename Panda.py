class Panda:

    def __init__(self, pos0, vit0, m, mass=1):
        """
        pos0 : np.array([qx,qy])
                Initial position
        vit0 = np.array([vx,vy])
                Initial velocity
        m = np.array([0,0,mz])
            Moment dipolaire
        """
        self.pos0 = pos0
        self.vit0 = vit0
        self.m = m
        self.mass=mass
        self.pos = pos0
        self.vit = vit0
        self.lastPos = None
        self.lastVit = None
        self.storePos = [pos0]
        self.storeVit = [vit0]

    def update_pos(self, newPos):
        """
        Store current position
        Set current position as last position
        Set new position as current position
        """
        self.lastPos = self.storePos[-1]
        self.pos = newPos
        self.storePos.append(self.pos)

    def update_vit(self, newVit):
        self.lastVit = self.storeVit[-1]
        self.vit = newVit
        self.storeVit.append(self.vit)

    def clear_panda(self):
        self.pos = self.pos0
        self.vit = self.vit0
        self.lastPos = None
        self.lastVit = None
        self.storePos = [self.pos0]
        self.storeVit = [self.vit0]
