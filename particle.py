class particle(object):
    def __init__(self, pos, vel, mass, charge, fixed = False):
        self.pos = pos
        self.prevPos = None
        self.vel = vel
        self.mass = mass
        self.charge = charge
        self.fixed = fixed
        self.acc = None
    
    def getMass(self):
        return self.mass
    def getCharge(self):
        return self.charge
    def getPos(self):
        return self.pos
    def getPrevPos(self):
        return self.prevPos
    def getVel(self):
        return self.vel
    def getAcc(self):
        return self.acc
    
    def setPos(self, pos):
        self.prevPos, self.pos = self.pos, pos
    def setVel(self, vel):
        self.vel = vel
    def setAcc(self, acc):
        self.acc = acc