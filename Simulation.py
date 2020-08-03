import random
import os
import numpy as np
import string
import pandas as pd
import MergingPolicies as mp
MAIN_DATA_DIRECTORY = '/Users/yishaiazabary/Desktop/My Career/OriAvinoamSimulation/Data/Simulation/BiggestFirst'
# MAIN_DATA_DIRECTORY = '/Users/yishaiazabary/Desktop/My Career/OriAvinoamSimulation/Data/Simulation/Random'


class Simulation:
    class CellSimulator:
        def __init__(self):
            self.nuclei_ctr = 1
            self.all_cells = []

        def add_cell(self, cell):
            self.nuclei_ctr += cell.nuclei_ctr
            self.all_cells.append(cell)

        def __repr__(self):
            return "NucleiCTR={}".format(self.nuclei_ctr)

    def __init__(self, initial_cell_number, merging_policy, simulation_name: str = ''):
        self.all_cells = [self.CellSimulator() for i in range(initial_cell_number)]
        self.mergin_policy = merging_policy
        if simulation_name is '':
            self.simulation_name = ''.join(random.choices(string.ascii_uppercase + string.digits, k=10))
        else:
            self.simulation_name = simulation_name

        self.simulation_df = pd.DataFrame(columns=['Frame#', 'CTRn=1','CTRn=2','CTRn=3','CTRn=4','CTR1total','CTR2total','CTR3total','CTR4total'])

    def update_df(self, frame_idx):
        # open a new row in the df
        self.simulation_df.loc[frame_idx, 'Frame#'] = frame_idx
        self.simulation_df.loc[
            frame_idx, ['CTRn=1', 'CTRn=2', 'CTRn=3', 'CTRn=4', 'CTR1total', 'CTR2total', 'CTR3total', 'CTR4total']] = 0
        # fill out the results
        for examined_cell in self.all_cells:
            number_of_cells_inside_examined_cell = len(examined_cell.all_cells) + 1
            number_of_nuclei_inside_examined_cell = examined_cell.nuclei_ctr
            if number_of_cells_inside_examined_cell == 1:
                col_idx = 1
                self.simulation_df.loc[frame_idx, 'CTRn={}'.format(col_idx)] += 1
                self.simulation_df.loc[frame_idx, 'CTR{}total'.format(col_idx)] += examined_cell.nuclei_ctr
            elif number_of_cells_inside_examined_cell == 2:
                col_idx = 2
                self.simulation_df.loc[frame_idx, 'CTRn={}'.format(col_idx)] += 1
                self.simulation_df.loc[frame_idx, 'CTR{}total'.format(col_idx)] += examined_cell.nuclei_ctr
            elif number_of_cells_inside_examined_cell == 3:
                col_idx = 3
                self.simulation_df.loc[frame_idx, 'CTRn={}'.format(col_idx)] += 1
                self.simulation_df.loc[frame_idx, 'CTR{}total'.format(col_idx)] += examined_cell.nuclei_ctr
            elif number_of_cells_inside_examined_cell >= 4:
                col_idx = 4
                self.simulation_df.loc[frame_idx, 'CTRn={}'.format(col_idx)] += 1
                self.simulation_df.loc[frame_idx, 'CTR{}total'.format(col_idx)] += examined_cell.nuclei_ctr

    def run_simulation_fusions_by_portion(self, frames_number: int = 16, portion_of_cell_to_merge: float = 0.25):
        for frame_idx in range(frames_number):
            ctr_merging_cells_in_frame = 0
            # go through a portion of the cells and see if they merge
            for cellIDX in range(int(len(self.all_cells) * portion_of_cell_to_merge)):
                # get target cell to merge into with merging policy
                idx_to_merge_into = np.random.randint(0, len(self.all_cells))
                indiv_cell = self.all_cells[idx_to_merge_into]
                target_cell = self.mergin_policy(self.all_cells)
                if target_cell is None:
                    continue
                ctr_merging_cells_in_frame += 1
                indiv_cell.add_cell(target_cell)
                self.all_cells.remove(target_cell)
            self.update_df(frame_idx)

    def run_simulation_fusions_by_fusions_list(self, n_fusions_per_frame: [int]):
        for frame_idx, n_fusions_in_frame in enumerate(n_fusions_per_frame):
            for fusion in range(n_fusions_in_frame):
                indiv_cell = self.all_cells[np.random.randint(0, len(self.all_cells))]
                target_cell = self.mergin_policy(self.all_cells)
                target_cell.add_cell(indiv_cell)
                self.all_cells.remove(indiv_cell)
            # ctrz_merging_cells_in_frame = 0
            # # go through a portion of the cells and see if they merge
            # for cellIDX in range(int(len(self.all_cells) * portion_of_cell_to_merge)):
            #     # get target cell to merge into with merging policy
            #     idx_to_merge_into = np.random.randint(0, len(self.all_cells))
            #     indiv_cell = self.all_cells[idx_to_merge_into]
            #     target_cell = self.mergin_policy(self.all_cells)
            #     if target_cell is None:
            #         continue
            #     ctr_merging_cells_in_frame += 1
            #     indiv_cell.add_cell(target_cell)
            #     self.all_cells.remove(target_cell)
            self.update_df(frame_idx)

    def convert_simulation_to_data_csv(self, path: str = ''):
        if path is '':
            csv_file_path = os.sep.join([MAIN_DATA_DIRECTORY, ''.join([self.simulation_name, '_', self.mergin_policy.__name__, '.csv'])])
        else:
            csv_file_path = path
        self.simulation_df.to_csv(csv_file_path)

#
# s = Simulation(initial_cell_number=200, merging_policy=mp.largest_first_merge)
# s.run_simulation_fusions_by_portion(frames_number=16)
# s.convert_simulation_to_data_csv()

# s = Simulation(initial_cell_number=200, merging_policy=mp.largest_first_merge)
# s.run_simulation_fusions_by_portion(frames_number=20)
# s.convert_simulation_to_data_csv()
