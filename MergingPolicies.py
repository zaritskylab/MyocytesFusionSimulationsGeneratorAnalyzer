import numpy as np
import pandas as pd
randomness_for_no_merge = -1


def random_merge(cell_list):
    no_merge_check = 0 if np.random.random_sample() <= randomness_for_no_merge else 1
    if no_merge_check == 1:
        return cell_list[np.random.randint(0, len(cell_list))]
    else:
        return None


def largest_first_merge(cell_list):
    no_merge_check = 0 if np.random.random_sample() <= randomness_for_no_merge else 1
    if no_merge_check == 1:
        max_cell = cell_list[0]

        for i in range(1, len(cell_list)):
            if cell_list[i].nuclei_ctr > max_cell.nuclei_ctr:
                max_cell = cell_list[i]

        return max_cell
    else:
        return None


def weighted_distribution_by_size(cell_list):
    dist_vec = np.zeros(len(cell_list))
    total_nuclei_ctr = 0
    max_nuclei_ctr = 0
    for cell in cell_list:
        if cell.nuclei_ctr > max_nuclei_ctr:
            max_nuclei_ctr = cell.nuclei_ctr
        total_nuclei_ctr += cell.nuclei_ctr
    for idx, cell in enumerate(cell_list):
        cell_score = float(cell.nuclei_ctr)/total_nuclei_ctr
        dist_vec[idx] = cell_score

    random_num = float(np.random.random_sample())

    return np.random.choice(cell_list, 1, p=dist_vec)[0]

    # for idx, dist in enumerate(dist_vec):
    #     if idx == 0:
    #         continue
    #     if random_num <= dist and random_num > dist[idx - 1]:
    #         return cell_list[idx]
    # return None