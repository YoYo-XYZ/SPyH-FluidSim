import numpy as np
from neighbor_search import _get_neighbor_cells
class Physattr:
    def __init__(self):
        # particle properties
        self.smoothing_radius = 0.1
        self.h = self.smoothing_radius
        self.kernelfunction = 'bspline' # currently only bspline
        self.rho = 1.0 #reference density
        self.mass = 1.0

        # governing physics
        self.gamma = 7.0
        self.viscosity = 1.0
        self.g = 981.0
        self.sound_speed = 34300.0

class Simattr(Physattr):
    def __init__(self):
        super().__init__()
        # particle spatial variables
        self.v_init = np.array([0,0])
        self.center_position = np.array([0,0])
        self.x_center = self.center_position[0]
        self.y_center = self.center_position[1]
        self.particle_amount_x = 5
        self.particle_amount_y = 5

        # physcial boundary variables
        self.collision_coefficient = 0.95
        self.bound = [2.5,2.5]

        # simulation config variables
        self.Dt = 0.1
        self.substep = 10
        self.time = 0
        self.g_acceleration = np.array([0,-self.g])
class Simulation(Simattr):
    def __init__(self):
        super().__init__()

    def initialize(self):
        self.particle_amount = self.particle_amount_x * self.particle_amount_y

        #position dynamics variables
        self.p_list = np.array([])
        self.v_list = np.array([self.v_init for i in range(0,self.particle_amount)])
        self.a_list = np.array([np.array([0,0]) for i in range(0,self.particle_amount)])
        self.rho_list = np.ones(self.particle_amount)
        self.g_list = np.array([self.g_acceleration for i in range(0,self.particle_amount)])

        self.realistic_parameter()
        self.setup_position()
        self.cell_size = 2*self.h
        self.cell_x = int(2*self.bound[0]/self.cell_size)

#-----------------Setup Initialize method-------------
    def realistic_parameter(self):
        self.mass = self.rho*4/3*np.pi*(2*self.h)**3
    def setup_position(self):
        self.x_gap = 1.99*self.h
        self.y_gap = 1.99*self.h

        self.p_list = np.array([self.center_position+np.array([self.x_gap*(x - ((self.particle_amount_x-1)/2)), self.y_gap*(y - ((self.particle_amount_y-1)/2)) ])
        for x in range(0,self.particle_amount_x)
        for y in range(0,self.particle_amount_y)])
        self.p_list_init = self.p_list

    def setup_random(self):
        self.p_list = np.random.rand(self.particle_amount,2)-np.array([[0.5,0.5]])

#-----------------Methods-------------------
    def distance_list_cal(self):
        self.index_map = _get_neighbor_cells(self.particle_amount, self.p_list, self.cell_size, self.cell_x)

        self.dvec_list = np.zeros([self.particle_amount,self.particle_amount,2])
        for i in range(0,self.particle_amount):
            for j in self.index_map[i]:
                if j>i:
                    self.dvec_list[i][j]=self.p_list[j]-self.p_list[i]
        #main variables
        self.dsize_list = np.linalg.norm(self.dvec_list,axis=2)

        #secondary variables
        self.r = self.dsize_list/self.h
        self.normvec = np.nan_to_num(self.dvec_list/self.dsize_list.reshape(self.r.shape + (1,)),nan=0)

    def kernel_cal(self):
        self.kernel_list = 3*(2*np.pi*self.h)*np.select([self.r==0,(self.r>0) & (self.r<1), (self.r>=1) & (self.r<2), self.r>=2], [0,2/3-self.r**2+1/2*self.r**3, (2-self.r)**3/6, 0])
        
    def grad_kernel_cal(self):
        x=np.select([self.r==0,(self.r>0) & (self.r<1), (self.r>=1) & (self.r<2), self.r>=2], [0,-2*self.r+3/2*self.r**2, -0.5*(2-self.r)**2, 0])
        self.grad_kernel_list = 3*(2*np.pi*self.h)*self.normvec*x.reshape(x.shape+(1,))

    def rho_div_cal(self):

        self.v_dif_list = np.zeros([self.particle_amount,self.particle_amount,2])
        for i in range(0,self.particle_amount):
            for j in self.index_map[i]:
                if j>i:
                    self.v_dif_list[i][j]=self.v_list[j]-self.v_list[i]

        self.v_dif_list = self.v_dif_list - np.transpose(self.v_dif_list,(1,0,2))

        self.alpha=self.v_dif_list/self.rho_list.reshape((1,)+self.rho_list.shape+(1,))

        self.x = np.einsum('ijk,ijk->ij', self.alpha, self.grad_kernel_list - np.transpose(self.grad_kernel_list,(1,0,2)))
        self.rho_div_list = self.mass*self.rho_list*np.sum(self.x,axis=1)

    def v_div_cal(self):
        #pressure force
        self.rho_av_list = np.zeros([self.particle_amount,self.particle_amount])
        for i in range(0,self.particle_amount):
            for j in self.index_map[i]:
                if j>i:
                    self.rho_av_list[i][j]=(self.sound_speed**2/self.gamma*(self.rho_list[i]**self.gamma+self.rho_list[j]**self.gamma-2))/(self.rho_list[i]*self.rho_list[j])


        x = self.rho_av_list.reshape(self.rho_av_list.shape+(1,))*self.grad_kernel_list
        x = x - np.transpose(x,axes=(1,0,2))#symmetize the value

        #viscosity force
        y = np.nan_to_num(np.linalg.norm(self.grad_kernel_list,axis = 2)/(self.dsize_list*self.rho_list.reshape((1,)+self.rho_list.shape)),nan=0)
        coif = 2*self.rho_list*self.viscosity
        yy=self.v_dif_list*y.reshape(y.shape+(1,))
        yy = yy - np.transpose(yy,(1,0,2))
        self.vis_force_list = coif.reshape(coif.shape+(1,))*np.sum(yy,axis=1)

        #net force
        self.v_div_list = self.mass*(np.sum(x,axis=1)+self.vis_force_list) + self.g_acceleration

#-----------------Time Integration Methods-------------------

    def integrate(self,dt):
        self.rho_list = self.rho_list + self.rho_div_list*dt

        #self.p_list = self.p_list + self.v_list*dt
        xxx = self.v_div_list*dt
        self.v_list = self.v_list + xxx
        self.p_list = self.p_list + self.v_list*dt + xxx*dt*0.5

    def update_position(self):

        dt=self.Dt/self.substep
        #loop substep

        for iii in range(0,self.substep):
            self.distance_list_cal()
            self.kernel_cal()
            self.grad_kernel_cal()
            self.rho_div_cal()
            self.v_div_cal()
            #updating actual values
            #print("a",self.grad_kernel_list)
            self.integrate(dt=dt)

            #checking box collision : x component
            abs_x = np.abs(self.p_list[:, 0])
            out_of_bounds = abs_x > self.bound[0]

            # Update velocity and position for out-of-bounds particles
            self.v_list[out_of_bounds, 0] *= -self.collision_coefficient
            self.p_list[out_of_bounds, 0] = np.copysign(self.bound[0], self.p_list[out_of_bounds, 0])

            #checking box collision : y component
            abs_y = np.abs(self.p_list[:, 1])
            out_of_bounds = abs_y > self.bound[1]

            # Update velocity and position for out-of-bounds particles
            self.v_list[out_of_bounds, 1] *= -self.collision_coefficient
            self.p_list[out_of_bounds, 1] = np.copysign(self.bound[1], self.p_list[out_of_bounds, 1])

        positions = self.p_list

        return positions
    

    def update_position_full(self,frame):
        position_superlist=[]
        velocitysize_superlist=[]
        dt=self.Dt/self.substep
        for i in range(0,frame):
            #loop substep
            for iii in range(0,self.substep):
                self.distance_list_cal()
                self.kernel_cal()
                self.grad_kernel_cal()
                self.rho_div_cal()
                self.v_div_cal()
                #updating actual values
                #print("a",self.grad_kernel_list)
                self.integrate(dt=dt)

                #checking box collision : x component
                abs_x = np.abs(self.p_list[:, 0])
                out_of_bounds = abs_x > self.bound[0]

                # Update velocity and position for out-of-bounds particles
                self.v_list[out_of_bounds, 0] *= -self.collision_coefficient
                self.p_list[out_of_bounds, 0] = np.copysign(self.bound[0], self.p_list[out_of_bounds, 0])

                #checking box collision : y component
                abs_y = np.abs(self.p_list[:, 1])
                out_of_bounds = abs_y > self.bound[1]

                # Update velocity and position for out-of-bounds particles
                self.v_list[out_of_bounds, 1] *= -self.collision_coefficient
                self.p_list[out_of_bounds, 1] = np.copysign(self.bound[1], self.p_list[out_of_bounds, 1])

            position_superlist.append(self.p_list)
            velocitysize_superlist.append(np.linalg.norm(self.v_list,axis=1))

        return position_superlist,velocitysize_superlist

