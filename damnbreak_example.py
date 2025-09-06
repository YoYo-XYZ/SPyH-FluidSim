import time
from Simulation import Simulation
from Illustration import illustrate

start_time = time.time()
sim0 = Simulation()
sim0.particle_amount_x=20
sim0.particle_amount_y=20
sim0.smoothing_radius = 0.1
sim0.particle_radius = 0.099
sim0.p_center_init = [0*-4.9,0*-2.5]
sim0.VISCOSITY = 0.01
sim0.Dt = 0.002
sim0.GRAVITY = 981
sim0.REST_PRESSURE = 1
sim0.substep = 40
sim0.frame = 200
sim0.bound = [4.5,4.5]
sim0.is_setup = False
sim0.initialize()
sim0.update_fulltime()

end_time = time.time()
print(f"time: {end_time-start_time}")

illustrate(sim0.position_superlist,sim0.velocitysize_superlist,sim0.bound,sim0.frame)