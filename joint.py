'''
Class defining the joint components mechanical properties and that computes MoS

'''

class Material():
    def __init__(self, E, nu):
        self.E = E
        self.nu = nu
        self.G = E / (2*(1+nu))


class Part(Material):
    def __init__(self, material, thickness, hole_diameter):
        self.material = material
        self.thickness = thickness
        self.hole_diameter = hole_diameter



class Bolt(Material):

    pi = 3.14159265359

    def __init__(self, system, material, thread_size, pitch, coeff_use):

        self.pitch = pitch
        self.h =  (3**(0.5))/2 * pitch #height of fundamental triangle
        #metric system, ext and int diameter in meters
        self.d_min_ext = (thread_size - 1.226869*pitch)        
        self.d_pitch_ext = (thread_size - 0.64951905*pitch)
        self.d_min_int = (thread_size - 1.08253175*pitch)
        self.d_pitch_int = (thread_size - 0.64951905*pitch)

        self.area_nom = self.pi*thread_size**2 / 4
        self.area_t = self.pi/4 * (thread_size -0.9382*pitch)**2
        self.area_min = self.pi/4 * self.d_min_ext**2


class Joint(Part, Bolt):
    def __init__(self, bolt, part_male, part_female, Type):
        self.bolt = bolt
        self.part_male = part_male
        self.part_female = part_female