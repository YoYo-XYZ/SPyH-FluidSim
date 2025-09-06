import numpy as np
from scipy.sparse import csr_array
from scipy.spatial import KDTree
import Attributes

def symmetrize_vector(array):
    """Return the antisymmetric part of a matrix."""
    return array - array.T
def symmetrize_scalar(array):
    """Return the symmetric part of a matrix."""
    return array + array.T

def bspline_kernel(R_data, kernel_constant):
    """Compute B-spline kernel and its gradient."""
    kernel_data = np.zeros_like(R_data)
    gradkernel_data = np.zeros_like(R_data)

    cond1 = (R_data > 0) & (R_data < 1)
    cond2 = (R_data >= 1) & (R_data < 2)

    kernel_data[cond1] = 2/3 - R_data[cond1]**2 + 1/2 * R_data[cond1]**3
    kernel_data[cond2] = (2 - R_data[cond2])**3 / 6

    gradkernel_data[cond1] = -2*R_data[cond1]+3/2*R_data[cond1]**2
    gradkernel_data[cond2] = -0.5*(2-R_data[cond2])**2

    return kernel_constant*kernel_data, kernel_constant*gradkernel_data


class Calculation(Attributes.SimulationAttributes):
    """Handles SPH calculations for a simulation step."""
    def __init__(self, sim_class, p_list, v_list, rho_list):
        super().__init__()
        self.__dict__.update(sim_class.__dict__)
        self.shape = (self.particle_amount, self.particle_amount)

        self.rho_list = rho_list
        self.p_list = p_list
        self.v_list = v_list

        self.x_list = self.p_list[:,0]
        self.y_list = self.p_list[:,1]

        self.vx_list = self.v_list[:,0]
        self.vy_list = self.v_list[:,1]

#------------------Primary Variables---------------------------------

    def _get_index_map(self):
        """Neighbor search using KDTree."""
        tree = KDTree(self.p_list)
        self.index_map = tree.query_ball_point(self.p_list, r=self.support_radius)

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
        
        kernel_data, gradkernel_data = bspline_kernel(R_ij_data, self.kernel_constant)

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
            self.a_bodyforce_list[:, 0] +
            self.mass * (symmetrize_vector(xxx) * self.rho_divide_list).sum(axis=1) * self.rho_divide_list
        )

        yyy = self.gradkernelsize_matrix * (
            self.P_ij_matrix * self.y_ij_norm_matrix -
            self.VISCOSITY * self.vy_ij_matrix / self.r_ij_matrix
        )
        vy_change = (
            self.a_bodyforce_list[:, 1] +
            self.mass * (symmetrize_vector(yyy) / self.rho_divide_list).sum(axis=1) * self.rho_divide_list
        )

        self.v_change_list = np.column_stack((vx_change, vy_change))

    def getallchange(self):
        """Compute all changes for a simulation step."""
        self._get_allprimaryval()
        self._get_rho_change()
        self._get_v_change()
        return self.rho_change_list, self.v_change_list

