import numpy as np
def count_particles(nparray):
    """Count the number of particles."""
    return len(nparray)

class ParticleInitialize():
    def __init__(self, gap):
        self.gap = gap
        self.p_list = np.empty((2,0))
        self.v_list = np.empty((2,0))

    def rectangle(self, width, length, center=[0,0], velocity=[0,0]):
        """Initialize particle positions in a rectangle."""
        x = np.arange(-length/2, length/2, self.gap) + center[0]
        y = np.arange(-width/2, width/2, self.gap) + center[1]
        X, Y = np.meshgrid(x, y)
        V_x, V_y = np.array([velocity[0] * np.ones_like(X), velocity[1] * np.ones_like(Y)])

        coords = np.vstack((X.ravel(), Y.ravel()))
        self.p_list = np.hstack((self.p_list, coords))
        
        self.v_list = np.hstack((self.v_list, np.vstack((V_x.ravel(), V_y.ravel()))))
        
    @property
    def particle_amount(self):
        return self.p_list.shape[1]
    @property
    def area_cover(self):
        return self.particle_amount * (self.gap**2)
    
if __name__ == "__main__":
    import matplotlib.pyplot as plt
    particles = ParticleInitialize(0.3)
    particles.rectangle(3, 3)
    particles.rectangle(3, 6, center=[3,3])
    plt.scatter(particles.p_list[0], particles.p_list[1], label=f'particles number : {particles.particle_amount}')
    plt.axis('equal')
    plt.legend()
    plt.show()