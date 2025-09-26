from .Calculation import Calculation
from .Attributes import SimAttributes, ParticleAttributes, PhysicsAttributes
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
class Simulation():
    """Main simulation class for SPH fluid."""
    def __init__(self, sim_config:SimAttributes, particle:ParticleAttributes):
        super().__init__()
        self.particle = particle
        self.is_setup = False
        self.__dict__.update(particle.__dict__)
        self.__dict__.update(sim_config.__dict__)

        self.timestepper = {'simple1':self._integrate1, 'simple2':self._integrate2, 'leapfrog':self._integrate_leapfrog}

    def _record_states(self):
        self.particle.p_list = self.p_list
        self.particle.v_list = self.v_list
        self.particle.rho_list = self.rho_list

    def _integrate1(self):
        calc = Calculation(self.particle, self.p_list, self.v_list, self.rho_list)
        rho_change_list, v_change_list = calc.getallchange()

        self.rho_list = self.rho_list + rho_change_list * self.dt
        dv = v_change_list * self.dt
        self.v_list = self.v_list + dv
        self.p_list = self.p_list + self.v_list * self.dt + dv * self.dt * 0.5

    def _integrate2(self):
        p_list_predict = self.p_list + self.v_list * self.dt
        calc = Calculation(self.particle, p_list_predict, self.v_list, self.rho_list)
        rho_change_list, v_change_list = calc.getallchange()

        self.rho_list = self.rho_list + rho_change_list * self.dt
        self.v_list = self.v_list + v_change_list * self.dt
        self.p_list = self.p_list + self.v_list * self.dt

    def _integrate_leapfrog(self):
        # Step 1: Calculate current acceleration and update velocity to a half-step
        calc = Calculation(self.particle, self.p_list, self.v_list, self.rho_list)
        rho_change_list_current, v_change_list_current = calc.getallchange()

        v_list_half_step = self.v_list + v_change_list_current * 0.5 * self.dt
        
        # Step 2: Update position
        self.p_list = self.p_list + v_list_half_step * self.dt

        # Step 3: Update density using the continuity equation
        # Note: The code already does this via rho_change_list, which is consistent with the continuity approach
        self.rho_list = self.rho_list + rho_change_list_current * self.dt
        
        # Step 4: Re-calculate new acceleration based on the new positions and densities
        calc_new = Calculation(self.particle, self.p_list, v_list_half_step, self.rho_list)
        _, v_change_list_new = calc_new.getallchange()
        
        # Step 5: Complete the velocity update for the full time step
        self.v_list = v_list_half_step + v_change_list_new * 0.5 * self.dt

    def _check_wall_collision(self):
        """Reflect particles at the simulation boundary."""
        abs_x = np.abs(self.p_list[0])
        out_of_bounds_x = abs_x > self.bound[0]
        self.v_list[0][out_of_bounds_x] *= -self.collision_coefficient
        self.p_list[0][out_of_bounds_x] = np.copysign(self.bound[0], self.p_list[0][out_of_bounds_x])

        abs_y = np.abs(self.p_list[1])
        out_of_bounds_y = abs_y > self.bound[1]
        self.v_list[1][out_of_bounds_y] *= -self.collision_coefficient
        self.p_list[1][out_of_bounds_y] = np.copysign(self.bound[1], self.p_list[1][out_of_bounds_y])

    def update_fulltime(self, step_method = 'simple2'):
        """Run the simulation for the specified number of frames."""
        frames = self.frames
        substep = self.substep

        self.step_method = step_method
        if self.is_setup:
            substep = 1
            frames = 1

        self.p_list_record = []
        self.v_list_record = []
        self.rho_list_record = []

        for frame_idx in range(frames):
            print(frame_idx)
            for _ in range(substep):
                self._integrate2()
                self.timestepper[step_method]()
                self._check_wall_collision()

            self.p_list_record.append(self.p_list)
            self.v_list_record.append(self.v_list)
            self.rho_list_record.append(self.rho_list)

        _ = np.array(self.v_list_record)
        self.vsize_list_record = np.sqrt(_[:,0]**2 +_[:,1]**2)

    def create_anim(self):
        """
        Visualize SPH simulation results as an animated scatter plot.
        """

        fig, ax = plt.subplots()
        scat = ax.scatter(
            self.p_list_record[0][0], self.p_list_record[0][1],
            c=self.vsize_list_record[0], cmap='rainbow', s=5,
                    alpha=1 , vmin=0, vmax=np.sqrt(986*1.5*self.bound[1])
        )
        plt.colorbar(scat, label='Speed')
        plt.xlabel("x (m)")
        plt.ylabel("y (m)")

        ax.set_xlim(-self.bound[0], self.bound[0])
        ax.set_ylim(-self.bound[1], self.bound[1])
        ax.set_aspect('equal', 'box')

        time_text = ax.text(0.01, 0.95, '', transform=ax.transAxes, fontsize=12, alpha = 1,
                        bbox=dict(facecolor='white', alpha=0.0, edgecolor='none'))
        substep_text = ax.text(
                0.02, 0.93, f'timestep {self.Dt/self.substep}, {self.kernel_type} kernel, {self.step_method}',
                transform=ax.transAxes, fontsize=5, alpha = 1,
                bbox=dict(facecolor='white', alpha=0.0, edgecolor='none'))
        viscosity_text = ax.text(
                0.02, 0.91, f'Viscosity{self.VISCOSITY}, Soundspeed: {self.SOUND_SPEED}',
                transform=ax.transAxes, fontsize=5, alpha = 1,
                bbox=dict(facecolor='white', alpha=0.0, edgecolor='none'))

        def update(frame):
            scat.set_offsets(self.p_list_record[frame].T)
            scat.set_array(self.vsize_list_record[frame].T)
            current_time = frame*self.Dt
            time_text.set_text(f'Time: {current_time:.4f} s')

        # Create the animation
        anim = animation.FuncAnimation(fig, update, interval=self.Dt*1000, frames=self.frames)
        
        plt.show()

        return anim

    def save_anim(self, anim, name):
        anim.save(f'{name}.gif',fps=5*1/self.Dt, writer='pillow')
        return anim



