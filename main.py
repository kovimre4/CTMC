import random
import math
from graphics import *
K = 8987551787.3681764
ELECTRON_MASS = 9.10938356*10**(-31)
NUCLEUS_MASS = 1.66053904*10**(-27)
UNIT_CHARGE = 1.6021766208*10**(-19)
R = 5.291772104242293*10**(-11)
V = 2187691.2630042396

class Vector(object):
    def __init__(self,x,y,z):
        self.components = [x,y,z]
        
    def __str__(self):
        return "{" + str(self.getX()) + "; " + str(self.getY()) + "; " + str(self.getZ()) + "}"
        
    def __repr__(self):
        return "{" + str(self.getX()) + "; " + str(self.getY()) + "; " + str(self.getZ()) + "}"
    
    def __add__(self, other):
        x = 0
        y = 0
        z = 0
        if type(other) == int or type(other) == float:
            x = self.getX() + other
            y = self.getY() + other
            z = self.getZ() + other
        elif type(other).__name__ == "Vector":
            x = self.getX() + other.getX()
            y = self.getY() + other.getY()
            z = self.getZ() + other.getZ()
        else:
            raise ValueError("Only scalars and other vectors can be added to vectors")
        return Vector(x,y,z)
    
    def __sub__(self, other):
        x = 0
        y = 0
        z = 0
        if type(other) == int or type(other) == float:
            x = self.getX() - other
            y = self.getY() - other
            z = self.getZ() - other
        elif type(other).__name__ == "Vector":
            x = self.getX() - other.getX()
            y = self.getY() - other.getY()
            z = self.getZ() - other.getZ()
        else:
            raise ValueError("Only scalars and other vectors can be subtracted from vectors")
        return Vector(x,y,z)
    
    def __mul__(self, n):
        x = 0
        y = 0
        z = 0
        if type(n) == int or type(n) == float:
            x = self.getX() * n
            y = self.getY() * n
            z = self.getZ() * n
        else:
            raise ValueError("Only scalars can multiply together with vectors")
        return Vector(x,y,z)
    
    def __rmul__(self, n):
        return self*n
    
    def __truediv__(self, other):
        x = 0
        y = 0
        z = 0
        if type(other) == int or type(other) == float:
            x = self.getX() / other
            y = self.getY() / other
            z = self.getZ() / other
        else:
            raise ValueError("Only scalars can divide vectors")
        return Vector(x,y,z)
    
    def getX(self):
        return self.components[0]
    
    def getY(self):
        return self.components[1]
    
    def getZ(self):
        return self.components[2]
    
    def magnitude(self):
        return (self.getX()**2+self.getY()**2+self.getZ()**2)**(1/2)

    def unit(self):
        return self/self.magnitude()

class Particle(object):
    def __init__(self, pos, vel, mass, charge, fixed = False):
        self.pos = pos
        self.vel = vel
        self.mass = mass
        self.charge = charge
        self.fixed = fixed
        self.acc = 0
        particles.append(self)
    
    def update(self):
        if not self.fixed:
            self.pos = self.getPos() + self.getVel()*dt + self.getAcc()/2*dt**2
            self.vel = self.getVel() + self.getAcc()*dt
    
    def getMass(self):
        return self.mass
    def getCharge(self):
        return self.charge
    def getPos(self):
        return self.pos
    def getVel(self):
        return self.vel
    def getAcc(self):
        return self.acc
    def calcAcc(self):
        netForce = Vector(0,0,0)
        for particle in particles:
            if particle != self:
                relPos = self.getPos()-particle.getPos()
                netForce = netForce + (K*self.getCharge()*particle.getCharge()*relPos)/(relPos.magnitude()**3)
        self.acc = netForce/self.getMass()

class H(object):
    def __init__(self, pos, vel):
        self.pos = pos
        self.vel = vel
        velDir = randomUnitVector()
        axisDir = Vector(velDir.getX(), (-velDir.getX()**2 - velDir.getZ()**2)/velDir.getY(), velDir.getZ()).unit()
        self.nucleus = Particle(-1*axisDir*R/(1+NUCLEUS_MASS/ELECTRON_MASS)+pos, -1*velDir*V/(1+NUCLEUS_MASS/ELECTRON_MASS)+vel, NUCLEUS_MASS, UNIT_CHARGE)
        self.electron = Particle(axisDir*R/(1+ELECTRON_MASS/NUCLEUS_MASS)+pos, velDir*V/(1+ELECTRON_MASS/NUCLEUS_MASS)+vel, ELECTRON_MASS, -UNIT_CHARGE)

    def getNucleus(self):
        return self.nucleus
    def getElectron(self):
        return self.electron
    def getCenterOfMass(self):
        return (self.getNucleus().getPos()*self.getNucleus().getMass()+self.getElectron().getPos()*self.getElectron().getMass())/(self.getNucleus().getMass()+self.getElectron().getMass())
    def getVel(self):
        return (self.getNucleus().getVel()*self.getNucleus().getMass()+self.getElectron().getVel()*self.getElectron().getMass())/(self.getNucleus().getMass()+self.getElectron().getMass())
    def getEnergy(self):
        E = 0
        E += 1/2*self.getNucleus().getMass()*self.getNucleus().getVel().magnitude()**2
        E += 1/2*self.getElectron().getMass()*self.getElectron().getVel().magnitude()**2
        E += K*self.getNucleus().getCharge()*self.getElectron().getCharge()/(self.getElectron().getPos()-self.getNucleus().getPos()).magnitude()
        return E

def randomUnitVector():
    theta = random.uniform(0, 2*math.pi)
    z = random.uniform(-1,1)
    return Vector((1-z**2)**(1/2)*math.cos(theta), (1-z**2)**(1/2)*math.sin(theta), z)
        
def display():
    for particle in particles:
        circle = Circle(Point(particle.getPos().getX()*10**12*4+300, particle.getPos().getY()*10**12*4+300), 50/(-particle.getPos().getZ()*10**11+10))
        circle.setFill(color_rgb(255,255,255))
        circle.draw(canvas)
    
def step():
    for particle in particles:
        particle.calcAcc()
    for particle in particles:
        particle.update()
    display()

particles = []
elapsed = 0
dt = 10**(-19)
helium1 = H(Vector(0,0,0), Vector(400000,0,0))
helium2 = H(Vector(3*R,0,0), Vector(-400000,0,0))
initialE = helium1.getEnergy()

canvas = GraphWin("My Window", 1200, 600)
canvas.setBackground(color_rgb(0,0,0))
for i in range(10000):
    step()
    elapsed+=dt
    print(i, "/10000")
print("Rate of energy loss:", (helium1.getEnergy()-initialE)/elapsed, "J/s")
canvas.getMouse()   

