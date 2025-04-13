import tkinter as tk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from illustration import Illustration
import numpy as np

class Siminterface:
    def __init__(self, root):
        self.root = root
        self.root.title('SPH Fluid Simulation')
        self.root.geometry('1240x980')

        self.main_sim = Illustration()

        self.scale_sizex = tk.Scale(root, from_=0, to=10, length=200, resolution=0.1 ,orient='horizontal', label='box size x')
        self.scale_sizex.pack()
        self.scale_sizex.set(self.main_sim.bound[0])

        self.scale_sizey = tk.Scale(root, from_=0, to=10, length=200, resolution=0.1 ,orient='horizontal', label='box size y')
        self.scale_sizey.pack()
        self.scale_sizey.set(self.main_sim.bound[1])

        self.label_initposx = tk.Label(root, text='Start Position X')
        self.label_initposx.pack()
        self.entry_initposx = tk.Entry(root)
        self.entry_initposx.pack()
        self.entry_initposx.insert(0, '0')

        self.label_initposy = tk.Label(root, text='Start Position Y')
        self.label_initposy.pack()
        self.entry_initposy = tk.Entry(root)
        self.entry_initposy.pack()
        self.entry_initposy.insert(0, '0')

        self.scale_particleamountx = tk.Scale(root, from_=0, to=30, length=200, resolution=1 ,orient='horizontal', label='particle numbers x')
        self.scale_particleamountx.pack()
        self.scale_particleamountx.set(self.main_sim.particle_amount_x)
        self.scale_particleamounty = tk.Scale(root, from_=0, to=30, length=200, resolution=1 ,orient='horizontal', label='particle numbers y')
        self.scale_particleamounty.pack()
        self.scale_particleamounty.set(self.main_sim.particle_amount_y)

        #PHYSICS
        self.entry_smoothing = tk.Entry(root)
        self.entry_smoothing.pack()
        self.scale_smoothing = tk.Scale(root, from_=0, to=1, length=200, resolution=0.01 ,orient='horizontal', label='Smoothing Radius')
        self.scale_smoothing.pack()
        self.scale_smoothing.set(self.main_sim.smoothing_radius) 

        self.label_viscosity = tk.Label(root, text='viscosity')
        self.label_viscosity.pack()
        self.entry_viscosity = tk.Entry(root)
        self.entry_viscosity.pack()
        self.entry_viscosity.insert(0, str(self.main_sim.viscosity))

        self.label_gravity = tk.Label(root, text='gravity')
        self.label_gravity.pack()
        self.entry_gravity = tk.Entry(root)
        self.entry_gravity.pack()
        self.entry_gravity.insert(0, str(self.main_sim.g_acceleration[1]))

        self.label_frames = tk.Label(root, text='Animation Frames')
        self.label_frames.pack()
        self.entry_frames = tk.Entry(root)
        self.entry_frames.pack()
        self.entry_frames.insert(0, '100')

        self.scale_camspeed = tk.Scale(root, from_=0.25, to=2, length=200, resolution=0.25 ,orient='horizontal', label='Camera Speedx')
        self.scale_camspeed.pack()
        self.scale_camspeed.set(self.main_sim.camspeed) 

        self.enter_button1 = tk.Button(root, text = 'Perform Simulation', command=self.showplot)
        self.enter_button1.pack()

        self.enter_button2 = tk.Button(root, text = 'Show Setup', command=self.showsetup)
        self.enter_button2.pack()

    def getval(self):
        try:
            self.smoothingradius = float(self.scale_smoothing.get())
        except:
            self.smoothingradius = float(self.entry_smoothing.get())
        self.particleamountx = int(self.scale_particleamountx.get())
        self.particleamounty = int(self.scale_particleamounty.get())
        self.viscosity = float(self.entry_viscosity.get())
        self.gravity = float(self.entry_gravity.get())
        self.boxsizex = float(self.scale_sizex.get())
        self.boxsizey = float(self.scale_sizey.get())
        self.frames = int(self.entry_frames.get())
        self.camspeed = float(self.scale_camspeed.get())
        self.initposx = float(self.entry_initposx.get())
        self.initposy = float(self.entry_initposy.get())

    def showsetup(self):
        self.getval()

        self.assignvaluetosim()
        self.main_sim.initialize()
        self.setupfig = self.main_sim.illustrate_setup()

        if hasattr(self, 'canvas'):
            self.canvas.get_tk_widget().destroy()

        self.canvas = FigureCanvasTkAgg(self.setupfig, master=self.root)
        self.canvas.get_tk_widget().pack()
        self.canvas.draw()

    def showplot(self):
        self.getval()

        self.assignvaluetosim()
        self.main_sim.initialize()
        print(self.main_sim.viscosity)
        self.main_sim.illustrate()

    def assignvaluetosim(self):
        self.main_sim.particle_amount_x = self.particleamountx
        self.main_sim.particle_amount_y = self.particleamounty
        self.main_sim.bound = [self.boxsizex, self.boxsizey]
        self.main_sim.smoothing_radius = self.smoothingradius
        self.main_sim.viscosity = self.viscosity
        self.main_sim.frames = self.frames
        self.main_sim.g_acceleration[1] = self.gravity
        self.main_sim.camspeed = self.camspeed
        self.main_sim.center_position = np.array([self.initposx,self.initposy])