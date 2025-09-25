import numpy as np

class KernelFunction:
    def __init__(self, h):
        self.kernel_func = {'bspline':self.bspline_func, 'gaussian':self.quadratic_func}
        self.kernel_const = {'bspline':15/(7*np.pi*h**2), 'gaussian':1/(np.pi*h**2)}

    def bspline_func(self, R_data):
        """Compute B-spline kernel and its gradient."""
        kernel_data = np.zeros_like(R_data)
        gradkernel_data = np.copy(kernel_data)

        cond1 = (R_data > 0) & (R_data < 1)
        cond2 = (R_data >= 1) & (R_data < 2)

        kernel_data[cond1] = 2/3 - R_data[cond1]**2 + 1/2 * R_data[cond1]**3
        kernel_data[cond2] = (2 - R_data[cond2])**3 / 6

        gradkernel_data[cond1] = -2*R_data[cond1]+3/2*R_data[cond1]**2
        gradkernel_data[cond2] = -0.5*(2-R_data[cond2])**2

        return self.kernel_const['bspline']*kernel_data, self.kernel_const['bspline']*gradkernel_data

    def quadratic_func(self, R_data):
        """Compute Cubic-spline kernel and its gradient."""
        kernel_data = np.zeros_like(R_data)
        gradkernel_data = np.copy(kernel_data)

        cond1 = (R_data > 0) & (R_data < 3)

        kernel_data[cond1] = self.kernel_const['gaussian'] * np.exp(-R_data[cond1]**2)

        gradkernel_data[cond1] = kernel_data*(-2*R_data[cond1])

        return kernel_data, gradkernel_data