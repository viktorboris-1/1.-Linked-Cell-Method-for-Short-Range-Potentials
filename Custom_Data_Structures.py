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



