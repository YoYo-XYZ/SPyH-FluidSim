import numpy as np

class Simulation:
    def __init__(self):
        #spatial space variables
        self.center_position = np.array([0,0])
        self.x_center = self.center_position[0]
        self.y_center = self.center_position[1]
        self.particle_amount_x = 5
        self.particle_amount_y = 5
        self.particle_amount = self.particle_amount_x*self.particle_amount_y
        self.bound = [2.5,2.5]

        #physical property
        self.mass = 1
        self.smoothing_radius = 0.1
        self.viscosity = 1
        self.gamma = 7

        #simulation parameter variables
        self.Dt = 0.1
        self.substep = 10
        self.collision_coefficient = 0.95
        self.pressure_multiplier = 10
        self.time = 0

        self.g_acceleration = np.array([0,-9.8])
        self.v_init = np.array([0,0])

    def initialize(self):
        self.particle_amount = self.particle_amount_x * self.particle_amount_y

        #display property
        self.r_display = (self.mass/(4*np.pi))**(1/3)

        #position dynamics variables
        self.p_list = np.array([])
        self.v_list = np.array([self.v_init for i in range(0,self.particle_amount)])
        self.a_list = np.array([np.array([0,0]) for i in range(0,self.particle_amount)])
        self.rho_list = np.ones(self.particle_amount)
        self.g_list = np.array([self.g_acceleration for i in range(0,self.particle_amount)])

        self.setup_position()

    def realistic_parameter(self):
        self.mass = 4/3*np.pi*(2*self.smoothing_radius)**3/self.particle_amount
        self.x_gap = 1.98*self.smoothing_radius
        self.y_gap = 1.98*self.smoothing_radius
#-----------------Setup-----------------
    #set up the initial position
    def setup_position(self):
        self.x_gap = 1.98*self.smoothing_radius
        self.y_gap = 1.98*self.smoothing_radius

        self.p_list = np.array([self.center_position+np.array([self.x_gap*(x - ((self.particle_amount_x-1)/2)), self.y_gap*(y - ((self.particle_amount_y-1)/2)) ])
        for x in range(0,self.particle_amount_x)
        for y in range(0,self.particle_amount_y)])
        self.p_list_init = self.p_list

    def setup_random(self):
        self.p_list = np.random.rand(self.particle_amount,2)-np.array([[0.5,0.5]])

    #type II boundary particlle
    def generate_boundary_particle(self):

        bound_p_list = np.array([np.array([self])])

#-----------------Linked List Algorithm-------------------
    def generate_key(self):
        self.key_list = np.round(self.p_list/(2*self.smoothing_radius)).astype(int)
        #print(self.key_list)

    def search_nearby_index(self, target_index):

        relative_positions = np.array([
            [0, 1], [0, -1], [1, 0], [-1, 0],
            [1, 1], [1, -1], [-1, 1], [-1, -1],
            [0, 0]
        ])
        
        target_key = self.key_list[target_index]
        nearby_keys = target_key + relative_positions
        
        # Use broadcasting to compare all keys at once
        mask = (self.key_list[:, None] == nearby_keys).all(axis=2)
        nearby_index_list = np.where(mask.any(axis=1))[0]


        return nearby_index_list.tolist()
    def nearby_index_list_cal(self):
        self.nearby_index_list = [self.search_nearby_index(j) for j in range(0,self.particle_amount)]

#-----------------Methods-------------------
    def distance_list_cal(self):
        
        self.generate_key()
        self.nearby_index_list_cal()


        self.dvec_list = np.zeros([self.particle_amount,self.particle_amount,2])
        for i in range(0,self.particle_amount):
            for j in self.nearby_index_list[i]:
                if j>i:
                    self.dvec_list[i][j]=self.p_list[j]-self.p_list[i]
        #main variables
        self.dsize_list = np.linalg.norm(self.dvec_list,axis=2)

        #secondary variables
        self.r = self.dsize_list/self.smoothing_radius
        self.normvec = np.nan_to_num(self.dvec_list/self.dsize_list.reshape(self.r.shape + (1,)),nan=0)

    def kernel_cal(self):
        self.kernel_list = 3*(2*np.pi*self.smoothing_radius)*np.select([self.r==0,(self.r>0) & (self.r<1), (self.r>=1) & (self.r<2), self.r>=2], [0,2/3-self.r**2+1/2*self.r**3, (2-self.r)**3/6, 0])
        
    def grad_kernel_cal(self):
        x=np.select([self.r==0,(self.r>0) & (self.r<1), (self.r>=1) & (self.r<2), self.r>=2], [0,-2*self.r+3/2*self.r**2, -0.5*(2-self.r)**2, 0])
        self.grad_kernel_list = 3*(2*np.pi*self.smoothing_radius)*self.normvec*x.reshape(x.shape+(1,))

    def rho_div_cal(self):

        self.v_dif_list = np.zeros([self.particle_amount,self.particle_amount,2])
        for i in range(0,self.particle_amount):
            for j in self.nearby_index_list[i]:
                if j>i:
                    self.v_dif_list[i][j]=self.v_list[j]-self.v_list[i]

        self.v_dif_list = self.v_dif_list - np.transpose(self.v_dif_list,(1,0,2))

        self.gamma=self.v_dif_list/self.rho_list.reshape((1,)+self.rho_list.shape+(1,))

        self.x = np.einsum('ijk,ijk->ij', self.gamma, self.grad_kernel_list - np.transpose(self.grad_kernel_list,(1,0,2)))
        self.rho_div_list = self.mass*self.rho_list*np.sum(self.x,axis=1)

    def v_div_cal(self):
        #pressure force

        self.rho_av_list = np.zeros([self.particle_amount,self.particle_amount])
        for i in range(0,self.particle_amount):
            for j in self.nearby_index_list[i]:
                if j>i:
                    self.rho_av_list[i][j]=(34300**2/1*(self.rho_list[i]**1+self.rho_list[j]**1-2))/(self.rho_list[i]*self.rho_list[j])


        x = self.rho_av_list.reshape(self.rho_av_list.shape+(1,))*self.grad_kernel_list
        x = x - np.transpose(x,axes=(1,0,2))#symmetize the value

        #viscosity force
        y = np.nan_to_num(np.linalg.norm(self.grad_kernel_list,axis = 2)/(self.dsize_list*self.rho_list.reshape((1,)+self.rho_list.shape)),nan=0)
        coif = self.mass*2*self.rho_list*self.viscosity
        yy=self.v_dif_list*y.reshape(y.shape+(1,))
        yy = yy - np.transpose(yy,(1,0,2))
        self.vis_force_list = coif.reshape(coif.shape+(1,))*np.sum(yy,axis=1)

        #net force
        self.v_div_list = self.mass*np.sum(x,axis=1) + self.g_acceleration + self.vis_force_list

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

