'''
Classes for mechanical joints

TODO

[] coord
    [] adopt a Nastran Card definition approach
        [?] inherit pyNastran content?
[] load vector transformation from coord_in to coor_out
    [] check pyNastran library
[] bolt
    [] focus on ECSS implementation 
        [] optionally US system
    [] compute bolt design params
    [] compute pre-load
        [] obtain seating torque
        [] bolt stress

[] joint
    [] compute mos
    [] compute joint stiffness
    [] thread shear 
        [] int and ext thread
    [] clamped parts stresses
        [] pull through
        [] bearing
        [] pin bearing

[] hardware size
    [] wrench size

'''

import numpy as np
import pandas as pd

frict_coeff_path = "_input/db_material/friction.csv"
frict_df = pd.read_csv(frict_coeff_path, skiprows=4).set_index(['material_A', 'material_B']) 

class Material():
    '''
    Defines material physical properties
    '''
    def __init__(self, reference, name, E, nu, rho):
        self.reference = reference  #reference to material source data
        self.name = name            #Alloy Name, e.g.Ti6Al4V
        self.E = E                  #Young Moduls, Pa
        self.nu = nu                #Poisson Ratio, -
        self.rho = rho              #Density, kg/m^3

        # G modulus for isotropic materials
        self.G = E / (2*(1+nu))     #Shear Modulus, Pa

class friction_coef():
    '''
    given materials pair, it returns min and max friction coefficient
    '''

    def __init__(self, material_a, material_b, df=frict_df):
        self.material_a = str(material_a)
        self.material_b = str(material_b)
        self.df = df

    def nus(self):
        'returns min, max friction coeff and data reference'
        nu_min = self.df.loc[self.material_a, self.material_b]['nu_min']
        nu_max = self.df.loc[self.material_a, self.material_b]['nu_max']
        reference = self.df.loc[self.material_a, self.material_b]['reference']
        nus = {'nu_min' : nu_min, 'nu_max' : nu_max, 'reference' : reference}
        return nus


    # def __repr__():
    #     return f"Friction coefficient data:\n{frict_df}"

class Coord:
    '''
    defines coordinate system as normalized array
    '''
    def __init__(self,cid, arr):
        self.cid = cid
        self.arr = arr / np.abs(arr).max(axis=0) #clmn normalized array



class Load():
    '''
    returns loads wrt. bolt coordinate system
    '''
    def __init__(self, load_in, coord_in, coord_out, pull_direction):
        self.load_in = load_in
        self.coord_in = coord_in
        self.coord_out = coord_out
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


class Joint(Part, Bolt):
    '''
    computes MoS at joint level
    '''
    def __init__(self, Type, bolt, washer, nut, clamp_part_male, clamp_part_female):
        self.bolt = bolt
        self.clamp_part_male = clamp_part_male      # clamped part with through hole, directly underneath bolt head / bolt washer
        self.clamp_part_female = clamp_part_female  # clamped part with threaded hole or at the bolt nut side
        self.type = type                            # bolt with nut or tapped joint

        def calc_mos_slippage(self):
            '''
            computes slippage MoS
            '''
            pass



# example, combining OOP with Pandas DataFrames (proof of concept)
# input data as dataframe used to create bolt, joints objects; MoS are obtained through objects methods

if __name__ == '__main__':

    AA7075 = Material('mat_db', 'AA7075',71e9,0.21,2100.)
    Ti6Al4V = Material('mat_db', 'Ti-6Al-4V',400e9,0.21,4600.) 
    nitronic = Material('mat_db','nitronic',201e9,0.21,4000.)
    inox = Material('mat_db', 'inox', 202e9,0.22, 3800.)

    coord_glob = Coord(1001,np.diag(np.array([1,1,1])))
    coord_local = Coord(1101, np.array([[0,-1,0],[1,0,0],[0,0,-1]]))

    df = pd.DataFrame([[nitronic,6e-3,1e-3,0.90, coord_glob, coord_glob,inox,AA7075],[Ti6Al4V,5e-3,1e-3,0.80, coord_glob,coord_local,inox,inox]], columns=['material','thread_size', 'pitch', 'coeff_usage', 'c_in', 'c_out','mat_washer', 'mat_thread'])
    df['bolt'] = df.apply(lambda X: Bolt('SI',X[0],X[1],X[2],X[3]), axis=1)
    df['bolt_name'] = df['bolt'].apply(lambda X: X.bolt_name)
    df['bolt_mat_name'] = df['material'].apply(lambda X: X.name)

    df['area_tensile'] = df['bolt'].apply(lambda X: X.area_t)
    df['cid_in'] = df['c_in'].apply(lambda X: X.cid)
    df['cid_out'] = df['c_out'].apply(lambda X: X.cid)
    
    # df['nus'] = df[['bolt','mat_washer',]].apply(lambda X: X.bolt_name)

    with pd.option_context('display.max_rows', None, 'display.max_columns', None):  #printing the entire array more options can be specified also
        print(df)