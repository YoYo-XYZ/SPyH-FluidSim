import numpy as np

class PhysicsAttributes():
    """Physical constants and fluid properties."""
    GRAVITY = 981.0
    VISCOSITY = 0.01
    WATER_DENSITY = 1.0
    SOUND_SPEED = 150000.0
    GAMMA = 7.0
    REST_PRESSURE = 0.0
    REST_DENSITY = WATER_DENSITY

    def __init__(self):
        self.PRESSURE_CONSTANT = self.REST_DENSITY * self.SOUND_SPEED ** 2 / self.GAMMA

class ParticleAttributes(PhysicsAttributes):
    """Particle and kernel configuration for SPH."""
    def __init__(self):
        super().__init__()
        # SPH
        self.smoothing_radius = 0.1
        self.support_radius = 2 * self.smoothing_radius
        self.particle_radius = self.smoothing_radius*0.99
        self.kerneltype = 'bspline'

        # number of particle
        self.particle_amount_x = 5
        self.particle_amount_y = 5
        self.particle_amount = self.particle_amount_x*self.particle_amount_y

        # value determining distribution of particle
        self.GAP_FACTOR = 1
        self.x_gap = 2 * self.GAP_FACTOR * self.particle_radius
        self.y_gap = 2 * self.GAP_FACTOR * self.particle_radius
        self.p_center_init = [0.0, 0.0]
        self.v_init = [0.0, 0.0]

    def _get_p_list(self):
        """Initialize particle positions in a grid."""
        self.p_list = np.array([
            np.array(self.p_center_init) + np.array([
                self.x_gap * (x - (self.particle_amount_x - 1) / 2),
                self.y_gap * (y - (self.particle_amount_y - 1) / 2)
            ])
            for x in range(self.particle_amount_x)
            for y in range(self.particle_amount_y)
        ])
    
    def _initialize_particle(self):
        """Initialize particle velocities and densities."""
        self.v_list = np.array([self.v_init for _ in range(self.particle_amount)])
        self.rho_list = np.full(self.particle_amount, self.WATER_DENSITY)
        self._get_p_list()


class SimAttributes(ParticleAttributes):
    """Simulation configuration and boundary settings."""
    def __init__(self):
        super().__init__()
        # Time stepping
        self.Dt = 0.1
        self.substep = 10
        self.time = 0.0
        self.dt = self.Dt / self.substep
        self.frame = 10

        # Boundary and collision
        self.bound = [5.0, 5.0]
        self.collision_coefficient = 0.999

class SimulationAttributes(SimAttributes):
    """Main simulation attribute class with initialization."""
    def __init__(self):
        super().__init__()

    def initialize(self):
        """Initialize all simulation parameters and arrays."""
        self.h = self.smoothing_radius
        self.support_radius = 2*self.smoothing_radius
        self.mass = self.WATER_DENSITY*np.pi*self.particle_radius**2
        self.particle_amount = self.particle_amount_x*self.particle_amount_y
        self.kernel_constant = 15/(7*np.pi*self.smoothing_radius**2)

        self.dt = self.Dt/self.substep

        #for computation
        self.blank_ij_list = np.zeros((self.particle_amount, self.particle_amount))
        self.blank_vector_ij_list = np.zeros((self.particle_amount, self.particle_amount, 2))

        self.a_bodyforce_list = np.array([[0,-self.GRAVITY] for _ in range(self.particle_amount)])

        self.x_gap = 2*self.GAP_FACTOR*self.particle_radius
        self.y_gap = 2*self.GAP_FACTOR*self.particle_radius
        self._initialize_particle()
    
