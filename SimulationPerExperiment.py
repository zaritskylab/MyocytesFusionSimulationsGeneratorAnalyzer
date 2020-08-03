import os
import string
import numpy as np
import pandas as pd
import Simulation as sml
from MergingPolicies import random_merge, largest_first_merge, weighted_distribution_by_size
MAIN_DATA_DIRECTORY_PATH = '/Users/yishaiazabary/PycharmProjects/MucsleCellFusionSimulation/Data'


class SimulationPerExperiment:
    def __init__(self, exp_prefix: str):
        self.exp_prefix = exp_prefix
        self.exp_path = os.sep.join([MAIN_DATA_DIRECTORY_PATH,exp_prefix[:-10], '{0}.csv'.format(exp_prefix)])
        self.exp_df = pd.read_csv(self.exp_path)
        self.max_number_of_cells = np.max(self.exp_df.loc[:, 'CTRn=1'].values)
        self.number_of_frames = len(self.exp_df.loc[:, 'CTRn=1'].values)

    def calc_experiment_fusions_in_frame(self):
        n_fusions_per_frame = [0]
        # cols = ['CTRn=1', 'CTRn=2', 'CTRn=3', 'CTRn=4']
        cols = ['CTRn=2', 'CTRn=3', 'CTRn=4']

        frame_idx = 0
        for frame_idx in range(len(self.exp_df.loc[0:, ['CTRn=1', 'CTRn=2', 'CTRn=3', 'CTRn=4']].values) - 1):
            fusions_in_frame = 0
            # frame_ctrs = {
            #     cols[0]: self.exp_df.loc[frame_idx, 'CTRn=1'],
            #     cols[1]: self.exp_df.loc[frame_idx, 'CTRn=2'],
            #     cols[2]: self.exp_df.loc[frame_idx, 'CTRn=3'],
            #     cols[3]: self.exp_df.loc[frame_idx, 'CTRn=4']
            # }
            # next_frame_ctrs = {
            #     cols[0]: self.exp_df.loc[frame_idx + 1, 'CTRn=1'],
            #     cols[1]: self.exp_df.loc[frame_idx + 1, 'CTRn=2'],
            #     cols[2]: self.exp_df.loc[frame_idx + 1, 'CTRn=3'],
            #     cols[3]: self.exp_df.loc[frame_idx + 1, 'CTRn=4']
            # }
            frame_ctrs = {
                cols[0]: self.exp_df.loc[frame_idx, 'CTRn=2'],
                cols[1]: self.exp_df.loc[frame_idx, 'CTRn=3'],
                cols[2]: self.exp_df.loc[frame_idx, 'CTRn=4']
            }
            next_frame_ctrs = {
                cols[0]: self.exp_df.loc[frame_idx + 1, 'CTRn=2'],
                cols[1]: self.exp_df.loc[frame_idx + 1, 'CTRn=3'],
                cols[2]: self.exp_df.loc[frame_idx + 1, 'CTRn=4']
            }
            for idx, ctr_col in enumerate(cols):
                num_of_nuclei = idx + 2
                f_n = (next_frame_ctrs[ctr_col] - frame_ctrs[ctr_col]) * (num_of_nuclei - 1)
                fusions_in_frame += 0 if f_n < 0 else f_n
            # for colIDX in range(1, len(cols)+1):
                # if initial_ctrs[cols[colIDX]] != next_frame[colIDX]:
                # f_n = int(initial_ctrs[cols[colIDX-1]] - next_frame[colIDX-1]) * (colIDX - 1)

                # f_n = initial_ctrs
                # fusions_in_frame += 0 if f_n < 0 else f_n
                # initial_ctrs[cols[colIDX]] = next_frame[colIDX]
            n_fusions_per_frame.append(fusions_in_frame)
        return n_fusions_per_frame


n_replicates = 1000
for exp_prefix in ['200203_S17_processed']:#, '200203_S19_processed', '200203_S22_processed', '200203_S24_processed', '200604_S06_processed','200604_S09_processed']:

# exp_prefix = '200203_S17_processed'
    spe = SimulationPerExperiment(exp_prefix)
    n_fusions_per_frame = spe.calc_experiment_fusions_in_frame()
    max_number_of_cells = spe.max_number_of_cells
    number_of_frames = spe. number_of_frames
    print('Number of Fusions: {0}\n Number of initial cells: {1}'.format(n_fusions_per_frame, max_number_of_cells))
    for i in range(n_replicates):
        print('{0}: replica No. {1}'.format(exp_prefix, i))
        # random
        simulation = sml.Simulation(initial_cell_number=max_number_of_cells,
                                    merging_policy=random_merge,
                                    simulation_name='{}_simulation'.format(exp_prefix))
        simulation.run_simulation_fusions_by_fusions_list(n_fusions_per_frame=n_fusions_per_frame)
        simulation.convert_simulation_to_data_csv(path=os.sep.join([MAIN_DATA_DIRECTORY_PATH, exp_prefix[:-10], '1000_randoms', '{0}_{1}_{2}.csv'.format(exp_prefix, random_merge.__name__, i)]))
        # weighted
        simulation = sml.Simulation(initial_cell_number=max_number_of_cells,
                                    merging_policy=weighted_distribution_by_size,
                                    simulation_name='{}_simulation'.format(exp_prefix))
        simulation.run_simulation_fusions_by_fusions_list(n_fusions_per_frame=n_fusions_per_frame)
        simulation.convert_simulation_to_data_csv(path=os.sep.join([MAIN_DATA_DIRECTORY_PATH, exp_prefix[:-10], '1000_weighted', '{0}_{1}_{2}.csv'.format(exp_prefix, weighted_distribution_by_size.__name__, i)]))
