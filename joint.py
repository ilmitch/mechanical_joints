'''
Classes for mechanical joints

TODO (feel free to help out!)

[] material, given an input csv file istantiate all material automatically (is it possible??)
[] coord
    [] adopt a Nastran Card definition approach
        [?] inherit pyNastran content?
[] load vector transformation from coord_in to coor_out
    [] check pyNastran library

[-] bolt
    [] compute bolt design params
    [] compute pre-load
        [] obtain seating torque
        [] bolt stress
    [] focus on ECSS implementation 
        [] optionally US system

[-] joint
    [-] compute mos
        [-] slippage
        [] others
    [] compute joint stiffness
    [] thread shear 
        [] int and ext thread
    [] clamped parts stresses
        [] pull through
        [] bearing
        [] pin bearing

[] hardware size
    [] wrench size

[-] implement __repr__ for all classes

'''

import numpy as np
import pandas as pd

frict_coeff_path = "_input/db_material/friction.csv"
frict_df = pd.read_csv(frict_coeff_path, skiprows=4).set_index(['material_A', 'material_B']) 

class Analysis():
    '''
    Analysis Params as Class Variables
    '''
    km = 1.20            # mathematical model load factor
    kq = 1.25            # test qualification load factor
    kp = 1.00            # project or additional load factor
    kld = 1.20           # local design safety factor, typical for fasteners
    sf_yld = 1.10        # yield safety factor
    sf_ult = 1.25        # ultimate safety factor
    compliance = {'min' : 0.20, 'max' : 0.40}    # clamped parts compliance factor

    level_limit = km * kq * kp * kld    #factor from limit level to design level
    level_qual = kp * kld               #factor from qualidication level to design level
    level_design = 1.0                  #factor from design level to design level, 1.0
    
    # mechanical engineering theory handbook 
    handbook = {'title' : "Space Engineering, Threaded Fasteners Handbook",
                'release_date' : "16 April 2010",
                'doc_number' : 'ECSS-E-HB-32-23A'}

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
        self.sig_yld = sig_yld      #Yield Strength, Pa
        self.sig_ult = sig_ult      #Ultimate Strength, Pa

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
        self.load_in = load_in  # input load, list
        self.pull_direction = pull_direction    # bolt pull direction, string, 'X' or 'Y' or 'Z'
        
        # expressing input load as bolt pull / shear load
        if self.pull_direction == 'X':
            self.pull = load_in[0]
            self.shear = (load_in[1]**2 + load_in[2]**2)**0.5 
        elif self.pull_direction == 'Y':
            self.pull = load_in[1]
            self.shear = (load_in[0]**2 + load_in[2]**2)**0.5 
        elif self.pull_direction == 'Z':
            self.pull = load_in[2]
            self.shear = (load_in[0]**2 + load_in[1]**2)**0.5 
        else:
            print('provided pull_direction is wrong!')
            raise Exception
        


class Part(Material):
    def __init__(self, material, thickness):
        self.material = material
        self.thickness = thickness


class Bolt(Analysis):
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
    '''
    Returns: Washer object
    '''
    def __init__(self, material):
        self.material = material

class Nut():
    '''
    Returns: Nut object
    '''
    def __init__(self, material):
        self.material = material

    def __repr__(self):
        return f"Name:\t{self.material.name}\n  E-moduls:\t{self.material.E}[Pa]\n  Poisson Ratio:\t{self.material.nu}[-]\n  Density:\t{self.material.rho}[kg/m^3]"

class Joint(Analysis):
    '''
    Returns: Washer object, computes MoS at joint level
    '''
    def __init__(self, type, bolt, washer, nut, clamp_part_male, clamp_part_female, load):
        self.type = type                            # bolt with nut or tapped joint, string: 'nut' or 'tap'
        self.bolt = bolt                            # Bolt object
        self.washer = washer                        # Bolt Washer Object
        self.nut = nut                              # Bolt Nut Object
        self.clamp_part_male = clamp_part_male      # clamped part with through hole, directly underneath bolt head / bolt washer
        self.clamp_part_female = clamp_part_female  # clamped part with threaded hole or at the bolt nut side
        self.load = load                            # Load object

        # retrieving part pairs friction values though Friction Class
        self.parts_frict = Friction(clamp_part_male.material, clamp_part_female.material).nus()


    def calc_mos_slippage_ult(self):
        '''
        computes ultimate slippage MoS, see Handbook, page §9.2.3 Friction Grip Strength Analysis, page 137
        '''
        mos_slip = ((self.bolt.pt_min - (1-self.compliance['min']) * self.load.pull) * self.parts_frict['nu_min']) / (self.load.shear * self.sf_ult)
        return mos_slip

    def __repr__(self):
        #returns joint data and mos summary
        return f"joint type: {self.type}\n bolt type: {self.bolt.name}\n bolt washer material: {self.washer.material.name}\n"


# example, combining OOP with Pandas DataFrames (proof of concept)
# input data as dataframe used to create bolt, joints objects; MoS are obtained through objects methods

if __name__ == '__main__':

    # dummy materials
    AA7075 = Material('mat_db', 'AA7075',71e9,0.21,2100.,300e6,400e6)
    Ti6Al4V = Material('mat_db', 'Ti-6Al-4V',400e9,0.21,4600.,800e6,900e6) 
    nitronic = Material('mat_db','nitronic',201e9,0.21,4000.,600e6,700e6)
    inox = Material('mat_db', 'inox', 202e9,0.22, 3800.,500e6,600e6) 


    # coords
    coord_glob = Coord(1001,np.diag(np.array([1,1,1])))
    coord_local = Coord(1101, np.array([[0,-1,0],[1,0,0],[0,0,-1]]))

    # Proof of Concept, input data through DataFrame:
    # this input data is used to instantiate differnt object directly in the DataFrame
    df = pd.DataFrame([[nitronic,6e-3,1e-3,0.90, coord_glob, coord_glob,inox,AA7075],[Ti6Al4V,5e-3,1e-3,0.80, coord_glob,coord_local,inox,inox]], columns=['material','thread_size', 'pitch', 'coeff_usage', 'c_in', 'c_out','mat_washer', 'mat_thread'])
    df['bolt'] = df.apply(lambda X: Bolt('SI',X[0],X[1],X[2],X[3]), axis=1)
    df['bolt_name'] = df['bolt'].apply(lambda X: X.name)
    df['bolt_mat_name'] = df['material'].apply(lambda X: X.name)

    df['area_tensile'] = df['bolt'].apply(lambda X: X.area_t)
    df['cid_in'] = df['c_in'].apply(lambda X: X.cid)
    df['cid_out'] = df['c_out'].apply(lambda X: X.cid)
    
    with pd.option_context('display.max_rows', None, 'display.max_columns', None):  #printing the entire array more options can be specified also
        print(df)

    # istantiation example
    bolt_m6 = Bolt('SI', nitronic, 6e-3, 1e-3, 0.90)
    washer_inox = Washer(inox)
    nut_inox = Nut(inox)
    cbush_force = Load([1000,250,250])
    part_a = Part(AA7075,5e-3)
    part_b = Part(AA7075,10e-3)
    j1 = Joint('nut', bolt_m6, washer_inox, nut_inox,part_a, part_b, cbush_force)  
    print(j1)