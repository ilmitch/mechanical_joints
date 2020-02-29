'''
Classes for mechanical joints

'''

import numpy as np
import pandas as pd

class Material():
    def __init__(self, reference, name, E, nu, rho):
        self.reference = reference  #reference to material source data
        self.name = name            #Alloy Name, e.g.Ti6Al4V
        self.E = E                  #Young Moduls, Pa
        self.nu = nu                #Poisson Ratio, -
        self.rho = rho              #Density, kg/m^3

        # G modulus for isotropic materials
        self.G = E / (2*(1+nu))     #Shear Modulus, Pa


class Coord:
    '''
    defines coordinate system as normalized array
    '''
    def __init__(self,cid, arr):
        self.cid = cid
        self.arr = arr / np.abs(arr).max(axis=0) #clmn normalized array



class Load():
    '''
    returns loads wrt. bolt cid
    '''
    def __init__(self, load_in, coord_in, coord_out, pull_direction):
        self.load_in = load_in
        self.cid_in = cid_in
        self.cid_out = cid_out
        self.pull_direction = pull_direction


class Part(Material):
    def __init__(self, material, thickness):
        self.material = material
        self.thickness = thickness


class Bolt(Material):
    '''
    Bolts attributes, derived params and computation of Safety Margins
    '''
    pi = 3.14159265359

    def __init__(self, system, material, thread_size, pitch, coeff_use):

        self.system = system            #metric or US
        self.material = material        #material object
        self.thread_size = thread_size
        self.pitch = pitch
        self.h =  (3**(0.5))/2 * pitch  #height of fundamental triangle

        #metric system, ext and int diameter in meters
        # diameters
        ## min diameters
        self.d_min_ext = (thread_size - 1.226869*pitch)        
        self.d_min_int = (thread_size - 1.08253175*pitch)

        ## pitch diamters
        self.d_pitch_ext = (thread_size - 0.64951905*pitch)
        self.d_pitch_int = (thread_size - 0.64951905*pitch)

        # areas; nominal, tensile stress, minimal
        self.area_nom = self.pi/4 * (thread_size)**2
        self.area_t = self.pi/4 * (thread_size -0.9382*pitch)**2
        self.area_min = self.pi/4 * (self.d_min_ext)**2
        # bolt name
        self.bolt_name = f"M{thread_size*1e3}x{pitch*1e3}"

        def calc_mos(self):
            pass


class Joint(Part, Bolt):
    def __init__(self, bolt, part_male, part_female, Type):
        self.bolt = bolt
        self.part_male = part_male
        self.part_female = part_female


# example, combining OOP with Pandas DataFrames (proof of concept)
import pandas as pd


Ti6Al4V = Material('Material_DB', 'Ti-6Al-4V',71e9,0.21,4600) 
coord_glob = Coord(1001,np.diag(np.array([1,1,1])))
coord_local = Coord(1101, np.array([[0,-1,0],[1,0,0],[0,0,-1]]))

df = pd.DataFrame([[Ti6Al4V,6e-3,1e-3,0.90, coord_glob, coord_glob],[Ti6Al4V,5e-3,1e-3,0.80, coord_glob,coord_local]], columns=['material','thread_size', 'pitch', 'coeff_usage', 'c_in', 'c_out'])
df['material_name'] = df['material'].apply(lambda X: X.name)
df['bolt'] = df.apply(lambda X: Bolt('SI',X[0],X[1],X[2],X[3]), axis=1)
df['area_tensile'] = df['bolt'].apply(lambda X: X.area_t)
df['bolt_name'] = df['bolt'].apply(lambda X: X.bolt_name)
df['cid_in'] = df['c_in'].apply(lambda X: X.cid)
df['cid_out'] = df['c_out'].apply(lambda X: X.cid)
df['bolt_name'] = df['bolt'].apply(lambda X: X.bolt_name)

with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
    print(df)