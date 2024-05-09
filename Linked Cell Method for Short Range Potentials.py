import math
import turtle
from Custom_Data_Structures import Particle, ParticleList


class Simulation:


    def __init__(self, parameter_input_file):
        input_file = open(parameter_input_file, 'r')

        self.DIM = int(input_file.readline().split(':')[1][:-1:])

        num_cells = input_file.readline().split(":")[1][:-1:].split(',')
        self.num_cells = [None] * 2
        self.num_cells[0] = int(num_cells[0])
        self.num_cells[1] = int(num_cells[1])

        dimensions = input_file.readline().split(":")[1][:-1:].split(',')
        self.dimensions = [None] * 2 
        self.dimensions[0] = int(dimensions[0])
        self.dimensions[1] = int(dimensions[1])

        self.r_cut = float(input_file.readline().split(":")[1])
        self.delta_t = float(input_file.readline().split(":")[1])
        self.t_end = float(input_file.readline().split(':')[1])
        self.sigma = float(input_file.readline().split(':')[1])
        self.epsilon = float(input_file.readline().split(':')[1])
        self.t_start = float(input_file.readline().split(":")[1])
        
        self.grid = []
        for x in range(self.num_cells[0]):
            for y in range(self.num_cells[1]):
                self.grid.append(ParticleList())

        self.sim_history = []


    def index(self, multi_index):
        if len(multi_index) == 2:
            return multi_index[0] + self.num_cells[0] * multi_index[1]
        else:
            return multi_index[0] + self.num_cells[0] * (multi_index[1] + self.num_cells[1] * multi_index[2])


    def particle_block(self, bottom_left_corner, velocity, mass, x_len, y_len):
        offset = (2 ** (1/6)) * self.sigma
        for x in range(x_len):
            for y in range(y_len):
                position = [bottom_left_corner[0] + (x * offset), bottom_left_corner[1] + (y* offset)]
                multi_index = []
                for d in range(len(position)):
                    multi_index.append(math.floor((position[d]/self.dimensions[d]) * self.num_cells[d]))
                self.grid[self.index(multi_index)].insert_particle(Particle(mass, position, velocity, [0,0], 2))
        

    def display_system(self):
        turtle.setup(500,500)
        turtle.speed("fastest")
        turtle.tracer(0,0)

        simulation = turtle.Turtle()
        simulation.hideturtle()
        simulation.color("blue")
        simulation.penup()
        
        for grid in self.sim_history:
            turtle.reset()
            for particle_list in grid:
                current_particle_node = particle_list.head
                while current_particle_node != None:
                    current_particle = current_particle_node.particle
                    print(current_particle.position)
                    simulation.goto(current_particle.position[0], current_particle.position[1])
                    simulation.dot()
                    current_particle_node = current_particle_node.next_particle
            turtle.update()        
        turtle.exitonclick()
        

    def draw_domain(self):
        turtle.setup(500,500)
        turtle.speed("fastest")
        turtle.tracer(0,0)
        #window = turtle.Screen()
        simulation = turtle.Turtle()
        simulation.hideturtle()
        simulation.color("blue")
        simulation.penup()
    
        for particle_list in self.grid:
            current_particle_node = particle_list.head
            while current_particle_node != None:
                current_particle = current_particle_node.particle
                print(current_particle.position)
                simulation.goto(current_particle.position[0], current_particle.position[1])
                simulation.dot()
                current_particle_node = current_particle_node.next_particle
    
        turtle.exitonclick()
            
        
    def time_integration_basis(self):
        self.compute_force_LC()
        while self.t_start < self.t_end:
            self.sim_history.append(self.grid)
            self.t_start += self.delta_t
            self.compute_position_LC()
            self.compute_force_LC()
            self.compute_velocity_LC()
            

    def compute_force_LC(self):
        offsets = [-1, 0, 1]
        for x in range(self.num_cells[0]):
            for y in range(self.num_cells[1]):
                cell_multi_index = [x, y]
                current_particle_node = self.grid[self.index(cell_multi_index)].head
                while current_particle_node is not None:
                    current_particle_node.particle.force = [0] * self.DIM
                    for x_offset in offsets:
                        for y_offset in offsets: 
                            neighbor_multi_index = [x + x_offset, y + y_offset]
                            
                            second_particle_node = self.grid[self.index(neighbor_multi_index)].head
                            while second_particle_node is not None:
                                if current_particle_node != second_particle_node:
                                    r = 0
                                    for d in range(self.DIM):
                                        r += (current_particle_node.particle.position[d] - second_particle_node.particle.position[d]) ** 2
                                    if r <= self.r_cut:
                                        self.force(current_particle_node, second_particle_node)   
                                second_particle_node = second_particle_node.next_particle
                    current_particle_node = current_particle_node.next_particle
    

    def force(self, particle1, particle2):
        # Lennard Jones Potential
        r = 0
        for d in range(self.DIM):
            r += (particle1.particle.position[d] - particle2.particle.position[d]) ** 2
        s = (self.sigma ** 2) / r
        f = 24 * self.epsilon * s / r * (1 - 2 * s)
        for d in range(self.DIM):
            particle1.particle.force[d] += f * (particle2.particle.position[d] - particle1.particle.position[d])


    def compute_position_LC(self):
        for x in range(self.num_cells[0]):
            for y in range(self.num_cells[1]):
                cell_multi_index = [x, y]
                current_particle_node = self.grid[self.index(cell_multi_index)].head
                while current_particle_node is not None:
                    self.update_position(current_particle_node.particle)
                    current_particle_node = current_particle_node.next_particle
        self.move_particles_LC()


    def update_position(self, particle):
        a = self.delta_t * 0.5 / particle.mass
        for d in range(self.DIM):
            particle.position[d] += self.delta_t * (particle.velocity[d] + a * particle.force[d])
            particle.old_force[d] = particle.force[d]


    def compute_velocity_LC(self):
        for x in range(self.num_cells[0]):
            for y in range(self.num_cells[1]):
                cell_multi_index = [x, y]
                current_particle_node = self.grid[self.index(cell_multi_index)].head
                while current_particle_node is not None:
                    self.update_velocity(current_particle_node.particle)
                    current_particle_node = current_particle_node.next_particle


    def update_velocity(self, particle):
        a = self.delta_t * 0.5 / particle.mass
        for d in range(self.DIM):
            particle.velocity[d] += a * (particle.force[d] + particle.old_force[d])


    def move_particles_LC(self):
        next_index = [None] * self.DIM
        for x in range(self.num_cells[0]):
            for y in range(self.num_cells[1]):
                current_particle_node = self.grid[self.index([x,y])].head
                while current_particle_node is not None:

                    for d in range(self.DIM):
                        next_index[d] = math.floor((current_particle_node.particle.position[d] / self.dimensions[d]) * self.num_cells[d])
                    
                    if next_index[0] != x or next_index[1] != y:
                            self.grid[self.index([x,y])].delete_particle(current_particle_node)
                            self.grid[self.index(next_index)].insert_particle(current_particle_node.particle)

                    current_particle_node = current_particle_node.next_particle
        

    def particle_cell_distance(self, current_particle_node, cell_index):
        cell_center_position = [0] * self.DIM
        for d in range(len(cell_center_position)):
            cell_center_position[d] = 0.5 + cell_index[d] * (self.dimensions[d]/self.num_cells[d])
            if self.DIM == 2:
                delta_x = current_particle_node.particle.position[0] - cell_center_position[0]
                delta_y = current_particle_node.particle.position[1] - cell_center_position[1]
                return math.sqrt((delta_x ** 2) + (delta_y ** 2))
            else:
                delta_x = current_particle_node.particle.position[0] - cell_center_position[0]
                delta_y = current_particle_node.particle.position[1] - cell_center_position[1]
                delta_z = current_particle_node.particle.position[2] - cell_center_position[2]
                return math.sqrt((delta_x ** 2) + (delta_y ** 2) + (delta_z ** 2))


    def isInBounds(self, particle_node):
        particle = particle_node.particle
        position = particle.position
        particle_position = particle.position
        for d in range(self.DIM):
            if particle_position[d] < 0 or particle_position[d] > self.dimensions[d]:
                return False
            else:
                return True


if __name__ == "__main__":
    sim = Simulation('Parameters.txt')
    sim.particle_block([200,240],[0,-10],1,4,4) 
    sim.particle_block([204,200],[0,0],1,20,20)
    sim.time_integration_basis(0)
    sim.draw_domain()