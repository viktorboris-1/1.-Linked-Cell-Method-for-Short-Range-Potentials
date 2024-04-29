class Particle:
    def __init__(self, mass, position, velocity, force):
        self.mass = mass
        self.position = position
        self.velocity = velocity
        self.force = force
        self.old_force = [None] * DIM


class ParticleNode:
    def __init__(self, particle):
        self.particle = particle
        self.next_particle = None


class ParticleList:
    def __init__(self):
        self.head = None
        self.count = 0

    def insert_particle(self, new_particle):
        new_particle_node = ParticleNode(new_particle)
        if self.head is None:
            self.head = new_particle_node
            return
        else:
            new_particle_node.next_particle = self.head
            self.head = new_particle_node

    def delete_particle(self):
        """TODO"""

class Simulation:
    DIM = 2
    num_cells = [None] * DIM
    input_file = open('Simulation_Parameters', 'r')
    params = input_file.readlines()
    N = params[0]
    pnc = 1
    dimensions = params[1].split(',')
    r_cut = params[2]
    delta_t = params[3]
    t_end = params[4]
    def __init__(self, parameter_input_file):
        input_file = open(parameter_input_file, 'r')
        self.DIM = input_file.readline().split(":")[1]
        self.num_cells = input_file.readline().split(":")[1].split()
        self.DIM = input_file.readline().split(":")[1]
        self.DIM = input_file.readline().split(":")[1]
        self.DIM = input_file.readline().split(":")[1]

