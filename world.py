import numpy as np
import constants as C
import particle as P
import graphics as gr


class world(object):
    def __init__(self, dt):
        self.dt = dt
        self.timeElapsed = 0
        self.particles = []
        self.usingVerletIntegration = False
    
    def eulerStep(self):
        """
        Stats were obtained from an array ofa 1000 simulations
        Mean of relative change in internal energy: 4.1288
        Median of relative change in internal energy: 0.0829
        Mean simulation time 7.2999s
        """
        self.setAccelerations()
        for particle in self.particles:
            particle.setPos(particle.getPos() + particle.getVel()*self.dt + 1/2*particle.getAcc()*self.dt**2 )
            particle.setVel(particle.getVel() + particle.getAcc()*self.dt)
        self.timeElapsed += self.dt
                
    def leapfrogStep(self):
        """
        Also called Velocity Verlet Algorithm
        Stats were obtained from an array ofa 1000 simulations
        Mean of relative change in internal energy: 0.9955
        Median of relative change in internal energy: 0.0258
        Mean simulation time 13.2813s
        """
        self.setAccelerations()
        newPos = []
        accAtNewPos = []
        for particle in self.particles:
            newPos.append(particle.getPos() + particle.getVel()*self.dt + 1/2*particle.getAcc()*self.dt**2 )
            accAtNewPos.append(self.calculateAccelerationAtPos(particle, newPos[-1]))
        for index, particle in enumerate(self.particles):
            particle.setPos(newPos[index])
            particle.setVel(particle.getVel() + 1/2*(particle.getAcc()+accAtNewPos[index])*self.dt)
        self.timeElapsed += self.dt
        
    def verletStep(self):
        """
        Stats were obtained from an array ofa 1000 simulations
        Mean of relative change in internal energy: 1.6383
        Median of relative change in internal energy: 7.7263e-06
        Mean simulation time 7.777s
        """
        self.usingVerletIntegration = True
        if self.timeElapsed == 0:
            self.eulerStep()
        else:
            self.setAccelerations()
            for particle in self.particles:
                particle.setPos(2*particle.getPos() - particle.getPrevPos() + particle.getAcc()*self.dt**2)
        self.timeElapsed += self.dt
    
    def setAccelerations(self):
        for particle in self.particles:
            particle.setAcc(self.calculateAccelerationAtPos(particle, particle.getPos()))
    
    def calculateAccelerationAtPos(self, particle, pos):
        netForce = np.array([0,0,0])
        for other in self.particles:
            if other != particle:
                relativePos = pos - other.getPos()
                netForce = netForce + C.COULOMBS_CONSTANT * (particle.getCharge()*other.getCharge()*relativePos) / (np.linalg.norm(relativePos)**3)
        return netForce / particle.getMass()
    
    def internalEnergy(self):
        if self.usingVerletIntegration:
            for particle in self.particles:
                particle.setVel((2*particle.getPos() - 2*particle.getPrevPos() + particle.getAcc()*self.dt**2)/(2*self.dt))
        potE = 0
        for particle in self.particles:
            for other in self.particles:
                if particle != other:
                    potE += 1/2 * C.COULOMBS_CONSTANT*particle.getCharge()*other.getCharge()/np.linalg.norm(particle.getPos()-other.getPos())
        kinE = 0
        for particle in self.particles:
            kinE += 1/2 * particle.getMass() * np.linalg.norm(particle.getVel())**2
        return potE + kinE
    
    def setInitialState(self):
        self.initialEnergy = self.internalEnergy()
        self.timeElapsed = 0
    
    def createHidrogen(self, pos, vel):
        velDir = self.randomUnitVector()
        axisDir = np.array([velDir[0], (-velDir[0]**2 - velDir[2]**2)/velDir[1], velDir[2]])
        axisDir = axisDir/np.linalg.norm(axisDir)
        self.particles.append(P.particle(-1*axisDir*C.R/(1+C.NUCLEUS_MASS/C.ELECTRON_MASS)+pos, -1*velDir*C.V/(1+C.NUCLEUS_MASS/C.ELECTRON_MASS)+vel, C.NUCLEUS_MASS, C.UNIT_CHARGE))
        self.particles.append(P.particle(axisDir*C.R/(1+C.ELECTRON_MASS/C.NUCLEUS_MASS)+pos, velDir*C.V/(1+C.ELECTRON_MASS/C.NUCLEUS_MASS)+vel, C.ELECTRON_MASS, -C.UNIT_CHARGE))
        
    def randomUnitVector(self):
        theta = np.random.uniform(0, 2*np.pi)
        z = np.random.uniform(-1,1)
        return np.array([(1-z**2)**(1/2)*np.cos(theta), (1-z**2)**(1/2)*np.sin(theta), z])
    
    def createStage(self, color = [255,255,255]):
        self.stage = gr.GraphWin("My Window", 1200, 600)
        self.stage.setBackground(gr.color_rgb(0,0,0))
        self.color = color
        
    def show(self):
        for particle in self.particles:
            pos = particle.getPos()
            if pos[2]*10**11<9:
                circle = gr.Circle(gr.Point(pos[0]*10**12*4+300, pos[1]*10**12*4+300), 50/(-pos[2]*10**11+10))
                circle.setFill(gr.color_rgb(self.color[0], self.color[1], self.color[2]))
                circle.draw(self.stage)
                
    def getInitialEnergy(self):
        return self.initialEnergy
    def getTimeElapsed(self):
        return self.timeElapsed