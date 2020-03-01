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

class Analysis():
    '''
    Analysis Params
    '''
    def __init__(self, km=1.20, kq=1.25, kp=1.0, sf_yld=1.10, sf_ult=1.25, sf_lcd=1.20):
        self.km = km            # mathematical model load factor
        self.kq = kq            # test qualificaiton load factor
        self.kp = kp            # project or additional load factor
        self.sf_yld = sf_yld    # yield safety factor
        self.sf_ult = sf_ult    # ultimate safety factor
        self.sf_lcd = sf_lcd    # local design safety factor

class Material():
    '''
    Defines physical properties of a material object
    '''
    def __init__(self, reference, name, E, nu, rho, sig_yld, sig_ult):
        self.reference = reference  #reference to material source data
        self.name = name            #Alloy Name, e.g.Ti6Al4V
        self.E = E                  #Young Moduls, Pa
        self.nu = nu                #Poisson Ratio, -
        self.rho = rho              #Density, kg/m^3
        self.sig_yld = sig_yld
        self.sig_ult = sig_ult

        # G modulus for isotropic materials
        self.G = E / (2*(1+nu))     #Shear Modulus, Pa

class Friction():
    '''
    given materials pair, it returns min and max friction coefficient
    the frict_df material names has to match with the material object name
    '''

    def __init__(self, material_a, material_b, df=frict_df):
        self.material_a = material_a.name
        self.material_b = material_b.name
        self.df = df

    def nus(self):
        'returns min, max friction coeff and data reference'
        nu_min = self.df.loc[self.material_a, self.material_b]['nu_min']
        nu_max = self.df.loc[self.material_a, self.material_b]['nu_max']
        reference = self.df.loc[self.material_a, self.material_b]['reference']
        nus = {'nu_min' : nu_min, 'nu_max' : nu_max, 'reference' : reference}
        return nus


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
    # def __init__(self, load_in, coord_in, coord_out, pull_direction):
    #     self.load_in = load_in
    #     self.coord_in = coord_in
    #     self.coord_out = coord_out
    #     self.pull_direction = pull_direction
    def __init__(self, load_in, pull_direction = 'Z'):
        self.load_in = load_in
        self.pull_direction = pull_direction


class Part(Material):
    def __init__(self, material, thickness):
        self.material = material
        self.thickness = thickness


class Bolt():
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
        self.name = f"M{thread_size*1e3}x{pitch*1e3}"

        #bolt pre-tension [] to be implemented as a function of bolt design and coeff of use
        self.pt_min = 1000  
        self.pt_nom = 1500
        self.pt_max = 2000


class Washer():
    def __init__(self, material):
        self.material = material

class Nut():
    def __init__(self, material):
        self.material = material

    def __repr__(self):
        return f"Name:\t{self.material.name}\n  E-moduls:\t{self.material.E}[Pa]\n  Poisson Ratio:\t{self.material.nu}[-]\n  Density:\t{self.material.rho}[kg/m^3]"

class Joint(Analysis):
    '''
    computes MoS at joint level
    '''
    def __init__(self, type, bolt, washer, nut, clamp_part_male, clamp_part_female, load):
        self.type = type                            # bolt with nut or tapped joint
        self.bolt = bolt                            # bolt object
        self.clamp_part_male = clamp_part_male      # clamped part with through hole, directly underneath bolt head / bolt washer
        self.clamp_part_female = clamp_part_female  # clamped part with threaded hole or at the bolt nut side
        self.parts_frict = Friction(clamp_part_male.material, clamp_part_female.material).nus()


        def calc_mos_slippage(self):
            '''
            computes slippage MoS
            '''
            pass



# example, combining OOP with Pandas DataFrames (proof of concept)
# input data as dataframe used to create bolt, joints objects; MoS are obtained through objects methods

if __name__ == '__main__':

    # materials
    AA7075 = Material('mat_db', 'AA7075',71e9,0.21,2100.,300e6,400e6)
    Ti6Al4V = Material('mat_db', 'Ti-6Al-4V',400e9,0.21,4600.,800e6,900e6) 
    nitronic = Material('mat_db','nitronic',201e9,0.21,4000.,600e6,700e6)
    inox = Material('mat_db', 'inox', 202e9,0.22, 3800.,500e6,600e6) 
    # coords
    coord_glob = Coord(1001,np.diag(np.array([1,1,1])))
    coord_local = Coord(1101, np.array([[0,-1,0],[1,0,0],[0,0,-1]]))

    # input data DataFrame
    df = pd.DataFrame([[nitronic,6e-3,1e-3,0.90, coord_glob, coord_glob,inox,AA7075],[Ti6Al4V,5e-3,1e-3,0.80, coord_glob,coord_local,inox,inox]], columns=['material','thread_size', 'pitch', 'coeff_usage', 'c_in', 'c_out','mat_washer', 'mat_thread'])
    df['bolt'] = df.apply(lambda X: Bolt('SI',X[0],X[1],X[2],X[3]), axis=1)
    df['bolt_name'] = df['bolt'].apply(lambda X: X.name)
    df['bolt_mat_name'] = df['material'].apply(lambda X: X.name)

    df['area_tensile'] = df['bolt'].apply(lambda X: X.area_t)
    df['cid_in'] = df['c_in'].apply(lambda X: X.cid)
    df['cid_out'] = df['c_out'].apply(lambda X: X.cid)
    
    with pd.option_context('display.max_rows', None, 'display.max_columns', None):  #printing the entire array more options can be specified also
        print(df)

    # general istantiation example
    bolt_m6 = Bolt('SI', nitronic, 6e-3, 1e-3, 0.90)
    washer_inox = Washer(inox)
    nut_inox = Nut(inox)
    cbush_force = Load([1000,250,250])
    part_a = Part(AA7075,5e-3)
    part_b = Part(AA7075,10e-3)
    j1 = Joint('nut', bolt_m6, washer_inox, nut_inox,part_a, part_b, cbush_force)  
