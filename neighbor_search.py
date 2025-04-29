import numpy as np

def _get_neighbor_cells(particle_amount, p_list, cell_size, cell_x):
    p_list = (p_list/cell_size).astype(np.int32)
    cell_map = p_list[:,0] + cell_x*p_list[:,1]

    index_map = {}

    for i in range(particle_amount):
        nearby_cell_list = []
        for dx in [-1, 0, 1]:
            for dy in [-1, 0, 1]:
                nearby_hash = cell_map[i] + (dx + dy * cell_x)
                nearby_cell_list.append(nearby_hash)

        index_list = []
        for j in range(particle_amount):
            if cell_map[j] in nearby_cell_list:
                index_list.append(j)

        index_map[i] = index_list

    return index_map

