import numpy as np
import constants as C
import world as W
import time


worlds = [W.world(10**(-19)) for i in range(1000)]
#world = W.world(10**(-19))
#world.createStage()

[world.createHidrogen(np.array([0,0,0]), np.array([400000,0,0])) for world in worlds]
[world.createHidrogen(np.array([3*C.R,0,0]), np.array([-400000,0,0])) for world in worlds]
#world.createHidrogen(np.array([0,0,0]), np.array([400000,0,0]))
#world.createHidrogen(np.array([3*C.R,0,0]), np.array([-400000,0,0]))

[world.setInitialState() for world in worlds]
#world.setInitialState()

start = time.time()
for i in range(10000):
    [world.eulerStep() for world in worlds]
#    world.leapfrogStep()
#    world.show()
    print(str(i)+"/10000")
end = time.time()

error = []
for world in worlds:
    val = (world.internalEnergy()-world.getInitialEnergy())/world.getInitialEnergy()
    error.append(val)
    print(val)


print("Mean of relative change in internal energy:", np.mean(np.abs(error)))
print("Median of relative change in internal energy:", np.median(np.abs(error)))
print("Mean simulation time", str((end-start)/len(error))+"s")