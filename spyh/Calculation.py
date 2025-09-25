import numpy as np
from scipy.sparse import csr_array
from scipy.spatial import KDTree
from .Attributes import PhysicsAttributes as ParticleAttributes

def symmetrize_vector(array):
    """Return the antisymmetric part of a matrix."""
    return array - array.T
def symmetrize_scalar(array):
    """Return the symmetric part of a matrix."""
    return array + array.T

class Calculation():
    """Handles SPH calculations for a simulation step."""
    def __init__(self, particle: ParticleAttributes, p_list=None, v_list=None, rho_list=None):
        self.__dict__.update(particle.__dict__)
        self.shape = (self.particle_amount, self.particle_amount)

        if p_list is not None:
            self.p_list = p_list
        if v_list is not None:
            self.v_list = v_list
        if rho_list is not None:
            self.rho_list = rho_list

        self.x_list = self.p_list[0]
        self.y_list = self.p_list[1]

        self.vx_list = self.v_list[0]
        self.vy_list = self.v_list[1]

#------------------Primary Variables---------------------------------

    def _get_index_map(self):
        """Neighbor search using KDTree."""
        tree = KDTree(self.p_list.T)
        self.index_map = tree.query_ball_point(self.p_list.T, r=self.support_radius)

    def _get_ij(self):
        """Get index pairs for neighbor interactions."""
        self._get_index_map()
        indices = np.array([
            (i, j)
            for i in range(self.particle_amount)
            for j in self.index_map[i] if j > i
        ])
        try:
            self.i_vals, self.j_vals = indices[:, 0], indices[:, 1]
        except Exception:
            self.i_vals, self.j_vals = [], []
        self.indices = (self.i_vals, self.j_vals)

    def _get_primaryval(self):
        """Compute primary SPH variables for all pairs."""
        vx_ij_data = self.vx_list[self.j_vals] - self.vx_list[self.i_vals]
        vy_ij_data = self.vy_list[self.j_vals] - self.vy_list[self.i_vals]

        x_ij_data = self.x_list[self.j_vals] - self.x_list[self.i_vals]
        y_ij_data = self.y_list[self.j_vals] - self.y_list[self.i_vals]
        r_ij_data = np.sqrt(x_ij_data**2 + y_ij_data**2)
        R_ij_data = r_ij_data/self.h
        
        kernel_data, gradkernel_data = self.kernel_func[self.kernel_type](R_ij_data)

        # r_ij_data[r_ij_data == 0] = 1e-9
        x_ij_norm_data = x_ij_data/r_ij_data
        y_ij_norm_data = y_ij_data/r_ij_data

        self.x_ij_norm_matrix = csr_array((x_ij_norm_data, self.indices), shape=self.shape)
        self.y_ij_norm_matrix = csr_array((y_ij_norm_data, self.indices), shape=self.shape)
        self.r_ij_matrix = csr_array((r_ij_data, self.indices), shape=self.shape)

        self.vx_ij_matrix = csr_array((vx_ij_data, self.indices), shape=self.shape)
        self.vy_ij_matrix = csr_array((vy_ij_data, self.indices), shape=self.shape)

        self.kernel_matrix = csr_array((kernel_data, self.indices), shape=self.shape)
        self.gradkernelsize_matrix = csr_array((gradkernel_data, self.indices), shape=self.shape)
    
    def _get_rho_divide(self):
        """Compute reciprocal of density."""
        self.rho_divide_list = 1/self.rho_list

    def _get_allprimaryval(self):
        self._get_ij()
        self._get_primaryval()
        self._get_rho_divide()

#---------------------------Governing Equation------------------
    def _get_P_ij(self):
        """Compute pressure for all pairs."""
        P_ij_data = (
            self.PRESSURE_CONSTANT *
            ((1 / self.REST_DENSITY) ** self.GAMMA *
             (self.rho_list[self.j_vals] ** self.GAMMA + self.rho_list[self.i_vals] ** self.GAMMA) - 2)
            + 2 * self.REST_PRESSURE
        )
        self.P_ij_matrix = csr_array((P_ij_data, self.indices), shape=self.shape)

    def _get_rho_change(self):
        """Compute density change for each particle."""
        x = self.rho_divide_list * symmetrize_scalar(
            (self.vx_ij_matrix * self.x_ij_norm_matrix + self.vy_ij_matrix * self.y_ij_norm_matrix)
            * self.gradkernelsize_matrix
        )
        self.rho_change_list = self.mass * self.rho_list * x.sum(axis=1)

    def _get_v_change(self):
        """Compute velocity change for each particle."""
        self._get_P_ij()

        xxx = self.gradkernelsize_matrix * (
            self.P_ij_matrix * self.x_ij_norm_matrix -
            self.VISCOSITY * self.vx_ij_matrix / self.r_ij_matrix
        )
        vx_change = (
            self.g_list[0] +
            self.mass * (symmetrize_vector(xxx) * self.rho_divide_list).sum(axis=1) * self.rho_divide_list
        )

        yyy = self.gradkernelsize_matrix * (
            self.P_ij_matrix * self.y_ij_norm_matrix -
            self.VISCOSITY * self.vy_ij_matrix / self.r_ij_matrix
        )
        vy_change = (
            self.g_list[1] +
            self.mass * (symmetrize_vector(yyy) * self.rho_divide_list).sum(axis=1) * self.rho_divide_list
        )

        self.v_change_list = np.array([vx_change, vy_change])

    def getallchange(self):
        """Compute all changes for a simulation step."""
        self._get_allprimaryval()
        self._get_rho_change()
        self._get_v_change()
        return self.rho_change_list, self.v_change_list

if __name__ == "__main__":
    from ParticleInit import ParticleInitialize
    points = ParticleInitialize(0.3)
    points.rectangle(5,5)
    particle = ParticleAttributes(0.3*2)
    particle._initialize_particle(points)

    import matplotlib.pyplot as plt
    plt.scatter(particle.p_list[0], particle.p_list[1], label=f'particles number : {particle.particle_amount}')
    plt.axis('equal')
    plt.legend()
    plt.show()

    cal = Calculation(particle)
    rho_change_list, v_change_list = cal.getallchange()
    print(cal.index_map)
    print(v_change_list)