import math
from Custom_Data_Structures import Particle, ParticleList


def main():
    DIM = 2
    simulation_dimensions = [None] * DIM


def compute_position_LC(grid, num_cells, dimensions, delta_t):
    for x in range(num_cells[0]):
        for y in range(num_cells[1]):
            cell_multi_index = [x, y]
            current_particle = grid[index(cell_multi_index, num_cells)]
            while current_particle is not None:
                update_position(current_particle.particle, delta_t)
                current_particle = current_particle.next_particle
    move_particles_LC(grid, num_cells, dimensions)


def compute_force_LC(grid, num_cells, cutoff_rad, DIM):
    offsets = [-1, 0, 1]
    for x in range(num_cells[0]):
        for y in range(num_cells[1]):
            cell_multi_index = [x, y]
            current_particle = grid[index(cell_multi_index, num_cells)]
            while current_particle is not None:
                particle_force = [0] * DIM
                for x_offset in offsets:
                    for y_offset in offsets:
                        neighbor_multi_index = [x + x_offset, y + y_offset]
                        # Treat kc[d] < 0 and kc[d] > nc[d] according to BCs
                        if particle_cell_distance(current_particle, neighbor_multi_index) <= cutoff_rad:
                            second_particle = grid[index(neighbor_multi_index, num_cells)]
                            while second_particle is not None:
                                if current_particle != second_particle:
                                    r = 0
                                    for d in range(DIM):
                                        r += (current_particle.particle.position - second_particle.particle.position) ** 2
                                        if r <= cutoff_rad:
                                            force(current_particle, second_particle)


def compute_velocity_LC(grid, num_cells, domain_dimensions, delta_t):
    for x in range(num_cells[0]):
        for y in range(num_cells[1]):
            cell_multi_index = [x, y]
            current_particle = grid[index(cell_multi_index, num_cells)]
            while current_particle is not None:
                update_velocity(current_particle.particle, delta_t)
                current_particle = current_particle.next_particle


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


def move_particles_LC(grid, num_cells, l, DIM):
    kc = [None] * DIM
    for x in range(num_cells[0]):
        for y in range(num_cells[1]):
            current_particle = grid(index([x,y],num_cells))
            while current_particle is not None:
                #treat BC for particle node
                for d in range(DIM):
                    kc[d] = math.floor(particle_list)



def particle_cell_distance(current_particle, neighbor_cell_index):
    """TODO"""


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


def output_results():
