import matplotlib.pyplot as plt
import matplotlib.animation as animation
#from new_function import Simulation


def illustrate(position_superlist, velocity_superlist, bound:list, frames):
    """
    Visualize SPH simulation results as an animated scatter plot.

    Args:
        position_superlist: List of particle position arrays for each frame.
        velocity_superlist: List of particle velocity magnitude arrays for each frame.
        bound: [x_bound, y_bound] for plot limits.
        frames: Number of frames to animate.
    """
    p_list = position_superlist[-1]
    v_list = velocity_superlist[-1]

    fig, ax = plt.subplots()
    scat = ax.scatter(
        p_list[:, 0], p_list[:, 1],
        c=v_list, cmap='rainbow', s=5, alpha=1 , vmin=0, vmax=10
    )
    plt.colorbar(scat, label='Speed')
    plt.xlabel("x (m)")
    plt.ylabel("y (m)")

    ax.set_xlim(-bound[0], bound[0])
    ax.set_ylim(-bound[1], bound[1])
    ax.set_aspect('equal', 'box')

    def update(frame):
        scat.set_offsets(position_superlist[frame])
        scat.set_array(velocity_superlist[frame])

    # Create the animation
    anim = animation.FuncAnimation(fig, update, interval=10, frames=frames)
    plt.show()
