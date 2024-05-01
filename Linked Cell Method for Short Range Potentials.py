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
        
        self.grid = []
        for x in range(self.num_cells[0]):
            for y in range(self.num_cells[1]):
                self.grid.append(ParticleList())

    
    def load_particles(self, particle_input_file):
        for x in range(self.num_cells[0]):
            for y in range(self.num_cells[1]):
                self.grid.append(ParticleList())

        input_file = open(particle_input_file, 'r')
        particles = input_file.readlines()

        for particle in particles:
            params = particle.split(',')
            mass = params[0]
            position = [params[1], params[2]]
            velocity = [params[3], params[4]]
            force = [params[5], params[6]]

            multi_index = []

            #Add 3D

            for d in range(len(position)):
                multi_index.append(math.floor((position[d]/self.dimensions[d]) * self.num_cells[d]))
            
            self.grid[self.index(multi_index, self.num_cells)] = ParticleList.insert_particle(Particle(mass, position, velocity, force, 2))
            

    def index(self, multi_index, num_cells):
        if len(multi_index) == 2:
            return multi_index[0] + num_cells[0] * multi_index[1]
        else:
            return multi_index[0] + num_cells[0] * (multi_index[1] + num_cells[1] * multi_index[2])
        
    def particle_block(self, sigma, top_left_corner, velocity, mass, x_len, y_len):
        offset = (2 ** (1/6)) * sigma
        multi_index = []

        for x in range(x_len):
            for y in range(y_len):
                position = [top_left_corner[0] + (x * offset), top_left_corner[1] + (y* offset)]
                for d in range(len(position)):
                    multi_index.append(math.floor((position[d]/self.dimensions[d]) * self.num_cells[d]))
            
                self.grid[self.index(multi_index, self.num_cells)].insert_particle(Particle(mass, position, velocity, [0,0], 2))

    def display_system(self):
        turtle.setup(400,400)
        turtle.speed("fastest")
        turtle.tracer(0,0)
        window = turtle.Screen()
        simulation = turtle.Turtle()
        simulation.hideturtle()
        simulation.color("blue")
        simulation.penup()
        
        for particle_list in self.grid:
            current_particle_node = particle_list.head
            while current_particle_node != None:
                current_particle = current_particle_node.particle
                simulation.goto(current_particle.position[0], current_particle.position[1])
                simulation.dot()
                current_particle_node = current_particle_node.next_particle
        turtle.update()
        
        window.exitonclick()

    def time_integration_basis(self, t, delta_t, t_end, grid, num_cells, domain_dimensions, cutoff_rad, DIM, sigma, epsilon):
        self.compute_force_LC(grid, num_cells, domain_dimensions, cutoff_rad, DIM, sigma, epsilon)
        while t < t_end:
            t += delta_t
            print("A")
            self.compute_position_LC(grid, num_cells, domain_dimensions, delta_t, DIM)
            print("B")
            self.compute_force_LC(grid, num_cells, domain_dimensions, cutoff_rad, DIM, sigma, epsilon)
            print("C")
            self.compute_velocity_LC(grid, num_cells, delta_t, DIM)
            print("D")
            #self.output_results()
    
    def compute_force_LC(self, grid, num_cells, domain_dimensions, cutoff_rad, DIM, sigma, epsilon):
        offsets = [-1, 0, 1]
        for x in range(num_cells[0]):
            for y in range(num_cells[1]):
                cell_multi_index = [x, y]
                current_particle_node = grid[self.index(cell_multi_index, num_cells)].head
                while current_particle_node is not None:
                   current_particle_node.particle.force = [0] * DIM
                   for x_offset in offsets:
                        for y_offset in offsets:
                            neighbor_multi_index = [x + x_offset, y + y_offset]
                            if x == num_cells[0] - 1 or y == num_cells[1] - 1: #BC
                              continue
                            if self.particle_cell_distance(current_particle_node, neighbor_multi_index, DIM, num_cells, domain_dimensions) <= cutoff_rad:
                                second_particle_node = grid[self.index(neighbor_multi_index, num_cells)].head
                                while second_particle_node is not None:
                                    if current_particle_node != second_particle_node:
                                        r = 0
                                        for d in range(DIM):
                                            r += (current_particle_node.particle.position[d] - second_particle_node.particle.position[d]) ** 2
                                            if r <= cutoff_rad:
                                                self.force(current_particle_node, second_particle_node, DIM, sigma, epsilon)   
                                    second_particle_node = second_particle_node.next_particle
                   current_particle_node = current_particle_node.next_particle
    

    def force(self, particle1, particle2, DIM, sigma, epsilon):
        # Lennard Jones Potential
        r = 0
        for d in range(DIM):
            r += (particle1.particle.position[d] - particle2.particle.position[d]) ** 2
        s = (sigma ** 2) / r
        f = 24 * epsilon * s / r * (1 - 2 * s)
        for d in range(DIM):
            particle1.particle.force[d] += f * (particle2.particle.position[d] - particle1.particle.position[d])


    def compute_position_LC(self, grid, num_cells, domain_dimensions, delta_t, DIM):
        for x in range(num_cells[0]):
            for y in range(num_cells[1]):
                cell_multi_index = [x, y]
                current_particle_node = grid[self.index(cell_multi_index, num_cells)].head
                while current_particle_node is not None:
                    self.update_position(current_particle_node.particle, delta_t, DIM)
                    current_particle_node = current_particle_node.next_particle
        self.move_particles_LC(grid, num_cells, domain_dimensions, DIM)


    def update_position(self, particle, delta_t, DIM):
        a = delta_t * 0.5 / particle.mass
        for d in range(DIM):
            particle.position[d] += delta_t * (particle.velocity[d] + a * particle.force[d])
            particle.old_force[d] = particle.force[d]


    def compute_velocity_LC(self, grid, num_cells, delta_t, DIM):
        for x in range(num_cells[0]):
            for y in range(num_cells[1]):
                cell_multi_index = [x, y]
                current_particle_node = grid[self.index(cell_multi_index, num_cells)].head
                while current_particle_node is not None:
                    self.update_velocity(current_particle_node.particle, delta_t, DIM)
                    current_particle_node = current_particle_node.next_particle


    def update_velocity(self, particle, delta_t, DIM):
        a = delta_t * 0.5 / particle.mass
        for d in range(DIM):
            particle.velocity[d] += a * (particle.force[d] + particle.old_force[d])


    def move_particles_LC(self, grid, num_cells, domain_dimensions, DIM):
        next_index = [None] * DIM
        for x in range(num_cells[0]):
            for y in range(num_cells[1]):
                current_particle_node = grid[self.index([x,y], num_cells)].head
                while current_particle_node is not None:
                    #if self.isInBounds(current_particle_node, self.dimensions, self.DIM) == False:
                        #current_particle_node = current_particle_node.next_particle
                        #grid[self.index([x,y], num_cells)].delete_particle(current_particle_node)
                        #continue

                    for d in range(DIM):
                        next_index[d] = math.floor((current_particle_node.particle.position[d] / domain_dimensions[d]) * num_cells[d])
                    
                    if next_index[0] != x or next_index[1] != y:
                        grid[self.index([x,y], num_cells)].delete_particle(current_particle_node)
                        grid[self.index(next_index, num_cells)].insert_particle(current_particle_node)
                    
                    current_particle_node = current_particle_node.next_particle


    def particle_cell_distance(self, current_particle_node, cell_index, DIM, num_cells, domain_dimensions):
        cell_center_position = [0] * DIM
        for d in range(len(cell_center_position)):
            cell_center_position[d] = 0.5 + cell_index[d] * (domain_dimensions[d]/num_cells[d])
            if DIM == 2:
                delta_x = current_particle_node.particle.position[0] - cell_center_position[0]
                delta_y = current_particle_node.particle.position[1] - cell_center_position[1]
                return math.sqrt((delta_x ** 2) + (delta_y ** 2))
            else:
                delta_x = current_particle_node.particle.position[0] - cell_center_position[0]
                delta_y = current_particle_node.particle.position[1] - cell_center_position[1]
                delta_z = current_particle_node.particle.position[2] - cell_center_position[2]
                return math.sqrt((delta_x ** 2) + (delta_y ** 2) + (delta_z ** 2))


    def isInBounds(self, particle_node, domain_dimensions, DIM):
        particle = particle_node.particle
        position = particle.position
        print(particle.position)
        particle_position = particle.position
        for d in range(DIM):
            if particle_position[d] < 0 or particle_position[d] > domain_dimensions[d]:
                return False
            else:
                return True


    def output_results():
        return

if __name__ == "__main__":
    #sim = Simulation(paramater_input_file)
    #sim.load_particles(particle_input_file)
    #time_integration_basis(0, sim.delta_t, sim.t_end, sim.grid, sim.num_cells, sim.dimensions, sim.r_cut, sim.DIM, sim.sigma, sim.epsilon)
    sim = Simulation('Parameters.txt')
    sim.particle_block(1,[75,100],[0,-10],1,40,40)
    sim.particle_block(1,[0,0],[0,-10],1,160,40)
    sim.display_system()
    sim.time_integration_basis(0, sim.delta_t, sim.t_end, sim.grid, sim.num_cells, sim.dimensions, sim.r_cut, sim.DIM, sim.sigma, sim.epsilon)
    sim.display_system()