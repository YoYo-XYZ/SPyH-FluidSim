from Calculation import Calculation
from Attributes import SimulationAttributes
import numpy as np

class Simulation(SimulationAttributes):
    """Main simulation class for SPH fluid."""
    def __init__(self):
        super().__init__()
        self.is_setup = False

    def _integrate(self):
        calc = Calculation(self, self.p_list, self.v_list, self.rho_list)
        rho_change_list, v_change_list = calc.getallchange()

        self.rho_list = self.rho_list + rho_change_list * self.dt
        uuu = v_change_list * self.dt
        self.v_list = self.v_list + uuu
        self.p_list = self.p_list + self.v_list * self.dt + uuu * self.dt * 0.5


    def _check_wall_collision(self):
        """Reflect particles at the simulation boundary."""
        abs_x = np.abs(self.p_list[:, 0])
        out_of_bounds_x = abs_x > self.bound[0]
        self.v_list[out_of_bounds_x, 0] *= -self.collision_coefficient
        self.p_list[out_of_bounds_x, 0] = np.copysign(self.bound[0], self.p_list[out_of_bounds_x, 0])

        abs_y = np.abs(self.p_list[:, 1])
        out_of_bounds_y = abs_y > self.bound[1]
        self.v_list[out_of_bounds_y, 1] *= -self.collision_coefficient
        self.p_list[out_of_bounds_y, 1] = np.copysign(self.bound[1], self.p_list[out_of_bounds_y, 1])

    def update_fulltime(self):
        """Run the simulation for the specified number of frames."""
        if self.is_setup:
            self.substep = 1
            self.frame = 1

        self.position_superlist = []
        self.velocitysize_superlist = []

        for frame_idx in range(self.frame):
            print(frame_idx)
            for _ in range(self.substep):
                # self._integrate_leapfrog()
                self._integrate()
                self._check_wall_collision()

            self.position_superlist.append(self.p_list.copy())
            self.velocitysize_superlist.append(np.linalg.norm(self.v_list, axis=1))


