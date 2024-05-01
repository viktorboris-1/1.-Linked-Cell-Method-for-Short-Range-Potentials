import math
import turtle

class Particle:
    def __init__(self, mass, position, velocity, force, DIM):
        self.mass = mass
        self.position = position
        self.velocity = velocity
        self.force = force
        self.old_force = [0] * DIM


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

    def delete_particle(self, particle_node_to_delete):
        current_node = self.head
        last_node = None
        if current_node == particle_node_to_delete:
            self.head = current_node.next_particle
            return
        
        while current_node != None and current_node != particle_node_to_delete:
            last_node = current_node
            current_node = current_node.next_particle

        if current_node is None:
            return
        else:
            last_node.next_particle = current_node.next_particle

