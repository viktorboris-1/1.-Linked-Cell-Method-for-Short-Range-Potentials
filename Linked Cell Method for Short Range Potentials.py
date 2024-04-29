import math
from Custom_Data_Structures import Particle, ParticleList


def main():
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

    for d in range(DIM)



def compute_position_LC(grid, num_cells, domain_dimensions, delta_t, DIM):
    for x in range(num_cells[0]):
        for y in range(num_cells[1]):
            cell_multi_index = [x, y]
            current_particle_node = grid[index(cell_multi_index, num_cells)]
            while current_particle_node is not None:
                update_position(current_particle_node.particle, delta_t, DIM)
                current_particle_node = current_particle_node.next_particle
    move_particles_LC(grid, num_cells, domain_dimensions, DIM)


def compute_force_LC(grid, num_cells, domain_dimensions, cutoff_rad, DIM, sigma, epsilon):
    offsets = [-1, 0, 1]
    for x in range(num_cells[0]):
        for y in range(num_cells[1]):
            cell_multi_index = [x, y]
            current_particle_node = grid[index(cell_multi_index, num_cells)]
            while current_particle_node is not None:
                current_particle_node.particle.force = [0] * DIM
                for x_offset in offsets:
                    for y_offset in offsets:
                        neighbor_multi_index = [x + x_offset, y + y_offset]
                        if x == len(num_cells[0]) - 1 or y == len(num_cells[1]) - 1: #BC
                            continue
                        if particle_cell_distance(current_particle_node, neighbor_multi_index, DIM, num_cells, domain_dimensions) <= cutoff_rad:
                            second_particle_node = grid[index(neighbor_multi_index, num_cells)]
                            while second_particle_node is not None:
                                if current_particle_node != second_particle_node:
                                    r = 0
                                    for d in range(DIM):
                                        r += (current_particle_node.particle.position - second_particle_node.particle.position) ** 2
                                        if r <= cutoff_rad:
                                            force(current_particle_node, second_particle_node, DIM, sigma, epsilon)


def compute_velocity_LC(grid, num_cells, domain_dimensions, delta_t, DIM):
    for x in range(num_cells[0]):
        for y in range(num_cells[1]):
            cell_multi_index = [x, y]
            current_particle_node = grid[index(cell_multi_index, num_cells)]
            while current_particle_node is not None:
                update_velocity(current_particle_node.particle, delta_t, DIM)
                current_particle_node = current_particle_node.next_particle


def update_position(particle, delta_t, DIM):
    a = delta_t * 0.5 / particle.mass
    for d in range(DIM):
        particle.position[d] += delta_t * (particle.velocity[d] + a * particle.force[d])
        particle.old_force[d] = particle.force[d]


def update_velocity(particle, delta_t, DIM):
    a = delta_t * 0.5 / particle.mass
    for d in range(DIM):
        particle.velocity[d] += a * (particle.force[d] + particle.old_force[d])


def index(multi_index, num_cells):
    if len(multi_index) == 2:
        return multi_index[0] + num_cells[0] * multi_index[1]
    else:
        return multi_index[0] + num_cells[0] * (multi_index[1] + num_cells[1] * multi_index[2])


def move_particles_LC(grid, num_cells, domain_dimensions, DIM):
    kc = [None] * DIM
    for x in range(num_cells[0]):
        for y in range(num_cells[1]):
            current_particle_node = grid(index([x,y],num_cells))
            while current_particle_node is not None:
                if isInBounds(current_particle_node) == False:

                for d in range(DIM):
                    kc[d] = math.floor(particle_list)



def particle_cell_distance(current_particle_node, cell_index, DIM, num_cells, domain_dimensions):
    cell_center_position = [0] * DIM
    for d in range(len(cell_center_position)):
        cell_center_position[d] = 0.5 + cell_index[d] * (domain_dimensions[d]/num_cells[d])
        if DIM == 2:
            delta_x = current_particle_node.particle.position[0] - cell_center_position[0]
            delta_y = current_particle_node.particle.position[1] - cell_center_position[1]
            return math.sqrt((delta_x ** 2) + (delta_y ** 2))
        else:
            return #TODO 3D



def force(particle1, particle2, DIM, sigma, epsilon):
    # Lennard Jones Potential
    r = 0
    for d in range(DIM):
        r += (particle1.particle.position - particle2.particle.position) ** 2
    s = (sigma ** 2) / r
    f = 24 * epsilon * s / r * (1 - 2 * s)
    for d in range(DIM):
        particle1 += f * (particle2.particle.position - particle1.particle.position)

def time_integration_basis(t, delta_t, t_end, grid, num_cells, domain_dimensions, cutoff_rad):
    compute_force_LC(grid, num_cells, cutoff_rad)
    while t < t_end:
        t += delta_t
        compute_position_LC(grid, num_cells, domain_dimensions, delta_t)
        compute_force_LC(grid, num_cells, cutoff_rad)
        compute_velocity_LC(grid, num_cells, domain_dimensions, delta_t)
        output_results()


def isInBounds(particle_node, domain_dimensions, DIM):
    particle_position = particle_node.particle.position
    for d in range(DIM):
        if particle_position[d] < 0 or particle_position > domain_dimensions[d]:
            return False
        else:
            return True


def output_results():
    """TODO"""