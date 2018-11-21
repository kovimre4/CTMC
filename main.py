import numpy as np
import constants as C
import world as W


#worlds = [W.world(10**(-19)) for i in range(10)]
world = W.world(10**(-19))
world.createStage()

#[world.createHidrogen(np.array([0,0,0]), np.array([400000,0,0])) for world in worlds]
#[world.createHidrogen(np.array([3*C.R,0,0]), np.array([-400000,0,0])) for world in worlds]
world.createHidrogen(np.array([0,0,0]), np.array([400000,0,0]))
world.createHidrogen(np.array([3*C.R,0,0]), np.array([-400000,0,0]))

#[world.setInitialState() for world in worlds]
world.setInitialState()

for i in range(10000):
    #[world.basicStep() for world in worlds]
    world.leapfrogStep()
    world.show()
    print(str(i)+"/10000")
    
    
print("Relative change in internal energy:", (world.internalEnergy()-world.getInitialEnergy())/world.getInitialEnergy())