import matplotlib.pyplot as plt
import matplotlib.animation as animation
from simulation_backend import Simulation
#from new_function import Simulation
import numpy as np


class Illustration(Simulation):
    def __init__(self):
        super().__init__()
        #Physics
        self.mass = 0.01
        self.viscosity = 10
        self.particle_amount_x = 10
        self.particle_amount_y = 10
        self.center_position = np.array([0,0])
        self.g_acceleration = np.array([0,-9.8])

        self.realistic_parameter()
        super().initialize()
        #Render
        self.camspeed = 1
        self.realtime_animation = False
        self.frames = 5
        self.substep = 20
        self.collision_coefficient = 0.9
        self.setup_position()

    def illustrate(self):
        self.Dt = 0.1*self.camspeed
        bound = self.bound
        fig, ax = plt.subplots()
        scat = ax.scatter(self.p_list[:, 0], self.p_list[:, 1], c=np.linalg.norm(self.v_list,axis=1), cmap='rainbow', s=567*self.smoothing_radius, alpha=1 , vmin=0, vmax=8)
        plt.colorbar(scat, label='Speed')
        fig.legend([f'viscosity = {self.viscosity}'])

        #animation computation
        ax.set_xlim(-bound[0], bound[0])
        ax.set_ylim(-bound[1], bound[1])
        ax.set_aspect('equal', 'box')

        if self.realtime_animation == True:
            def update(frame):
                position_list = self.update_position()
                scat.set_offsets(position_list)
                scat.set_array(np.linalg.norm(self.v_list,axis=1))

                return scat,
        else:
            position_superlist,velocitysize_superlist = self.update_position_full(self.frames)
            def update(frame):
                scat.set_offsets(position_superlist[frame])
                scat.set_array(velocitysize_superlist[frame])

                return scat,
    
        # Create the animation
        anim = animation.FuncAnimation(fig, update, interval=10, frames=self.frames)
        plt.show()

    def illustrate_setup(self):
        bound = self.bound
        fig, ax = plt.subplots()
        scat = ax.scatter(self.p_list[:, 0], self.p_list[:, 1], c=np.linalg.norm(self.v_list,axis=1), cmap='rainbow', s=567*self.smoothing_radius, alpha=1 , vmin=0, vmax=8)
        plt.colorbar(scat, label='Speed')
        fig.legend([f'viscosity = {self.viscosity}'])

        #animation computation
        ax.set_xlim(-bound[0], bound[0])
        ax.set_ylim(-bound[1], bound[1])
        ax.set_aspect('equal', 'box')

        scat.set_offsets(self.p_list_init)
        scat.set_array(np.linalg.norm(self.v_list,axis=1))

        return fig
