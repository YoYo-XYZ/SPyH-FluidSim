import numpy as np
from .ParticleInit import ParticleInitialize
from .KernelFunc import *
class PhysicsAttributes():
    def __init__(self, GRAVITY=981.0, VISCOSITY=30, WATER_DENSITY=1.0, SOUND_SPEED=1000.0, GAMMA=7.0, REST_PRESSURE=0.0):
        """Physical constants and fluid properties."""
        self.GRAVITY = GRAVITY
        self.VISCOSITY = VISCOSITY
        self.WATER_DENSITY = WATER_DENSITY
        self.SOUND_SPEED = SOUND_SPEED
        self.GAMMA = GAMMA
        self.REST_PRESSURE = REST_PRESSURE

    def calculate_constants(self):
        self.REST_DENSITY = self.WATER_DENSITY
        self.PRESSURE_CONSTANT = self.REST_DENSITY * self.SOUND_SPEED ** 2 / self.GAMMA
    
class ParticleAttributes(KernelFunction):
    """Particle and kernel configuration for SPH."""
    def __init__(self, support_radius, kernel_type = 'bspline'):
        # SPH
        self.support_radius = support_radius
        self.smoothing_radius = self.support_radius/2
        self.kernel_type = kernel_type
        KernelFunction.__init__(self, self.smoothing_radius)

    def _initialize_particle(self, particle:ParticleInitialize, physics:PhysicsAttributes):
        """Initialize particle velocities and densities."""
        physics.calculate_constants()
        self.__dict__.update(physics.__dict__)
        self.__dict__.update(particle.__dict__)
        self.particle_amount = particle.particle_amount
        self.rho_list = self.REST_DENSITY * np.ones(self.particle_amount)
        self.g_list = np.array([[0],[-self.GRAVITY]]) * np.ones_like(self.p_list)
        self.mass = self.REST_DENSITY*self.gap**2
        self.h = self.smoothing_radius

class SimAttributes():
    """Simulation configuration and boundary settings."""
    def __init__(self, bound:list, collision_coefficient, Dt, substep, totaltime):
        super().__init__()
        # Time stepping
        self.timestep_type = ["euler","verlet","leapfrog"]
        # Boundary and collision
        self.bound = bound
        self.collision_coefficient = collision_coefficient

        #visualization
        self.Dt = Dt
        self.substep = substep
        self.dt = self.Dt/self.substep

        self.totaltime = totaltime
        self.frames = int(self.totaltime//self.Dt)

if __name__ == "__main__":
    from .ParticleInit import ParticleInitialize
    points = ParticleInitialize(0.3)
    points.rectangle(5,5)
    particle = ParticleAttributes(0.3*1.3)
    particle._initialize_particle(points)

    import matplotlib.pyplot as plt
    plt.scatter(particle.p_list[0], particle.p_list[1], label=f'particles number : {particle.particle_amount}')
    plt.axis('equal')
    plt.legend()
    plt.show()