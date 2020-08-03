import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
MAIN_DATA_DIRECTORY_PATH = '/Users/yishaiazabary/PycharmProjects/MucsleCellFusionSimulation/Data'


class BootstrappingPerExperiment:
    def __init__(self, exp_prefix, type_of_simulation, col_to_calc_by, number_of_replications, frame_num):
        self.exp_prefix = exp_prefix
        self.type_of_simulation = type_of_simulation
        self.col_to_calc_by = col_to_calc_by
        type_of_simulation = ['{0}_{1}'.format(int(number_of_replications), type_of_simulation[0]),
                              '{0}_{1}'.format(int(number_of_replications), type_of_simulation[1])]
        self.path_to_simulations_directory = [os.sep.join([MAIN_DATA_DIRECTORY_PATH, exp_prefix, type_of_simulation[0]]),
                                              os.sep.join([MAIN_DATA_DIRECTORY_PATH, exp_prefix, type_of_simulation[1]])
                                              ]
        self.number_of_replications = number_of_replications
        self.frame_num = frame_num
        self.path_to_true_exp_file = os.sep.join([MAIN_DATA_DIRECTORY_PATH, exp_prefix, ''.join([exp_prefix, '_processed.csv'])])
        self.random_simulations_data_points = np.ndarray((self.number_of_replications, self.frame_num))
        self.weighted_simulations_data_points = np.ndarray((self.number_of_replications, self.frame_num))
        self.random_simulations_dist_averages = np.ndarray(self.frame_num)
        self.weighted_simulations_dist_averages = np.ndarray(self.frame_num)
        self.p_values_for_random_per_frame = np.zeros(frame_num)
        self.p_values_for_weighted_per_frame = np.zeros(frame_num)
        entire_exp_df = pd.read_csv(self.path_to_true_exp_file)
        # sum_of_totals = entire_exp_df.loc[:, 'CTRn=1'].values + \
        #                 entire_exp_df.loc[:, 'CTRn=2'].values + \
        #                 entire_exp_df.loc[:, 'CTRn=3'].values + \
        #                 entire_exp_df.loc[:, 'CTRn=4'].values
        self.initial_number_of_nuclei = np.max(entire_exp_df.loc[:, 'CTRn=1'].values)
        self.exp_data_points_normalized = entire_exp_df.loc[:, self.col_to_calc_by].values / self.initial_number_of_nuclei

    def read_simulations_col_data_points(self):
        random_files = (file for file in os.listdir(self.path_to_simulations_directory[0])
                 if os.path.isfile(os.path.join(self.path_to_simulations_directory[0], file)))
        weighted_files = (file for file in os.listdir(self.path_to_simulations_directory[1])
                        if os.path.isfile(os.path.join(self.path_to_simulations_directory[1], file)))
        idx = 0
        for filename in random_files:
            if filename.endswith('.csv'):
                path_to_single_file = os.sep.join([self.path_to_simulations_directory[0], filename])
                entire_df =  pd.read_csv(path_to_single_file)
                single_file_col_data =entire_df.loc[:, self.col_to_calc_by].values
                # sum_of_totals = entire_df.loc[:, 'CTRn=1'].values + \
                #                 entire_df.loc[:, 'CTRn=2'].values + \
                #                 entire_df.loc[:, 'CTRn=3'].values + \
                #                 entire_df.loc[:, 'CTRn=4'].values
                single_file_data_normalized = single_file_col_data / self.initial_number_of_nuclei
                self.random_simulations_data_points[idx] = single_file_data_normalized
                idx += 1
        idx = 0
        for filename in weighted_files:
            if filename.endswith('.csv'):
                path_to_single_file = os.sep.join([self.path_to_simulations_directory[1], filename])
                entire_df = pd.read_csv(path_to_single_file)
                single_file_col_data = entire_df.loc[:, self.col_to_calc_by].values
                # sum_of_totals = entire_df.loc[:, 'CTRn=1'].values + \
                #                 entire_df.loc[:, 'CTRn=2'].values + \
                #                 entire_df.loc[:, 'CTRn=3'].values + \
                #                 entire_df.loc[:, 'CTRn=4'].values
                # single_file_data_normalized = single_file_col_data / sum_of_totals
                self.weighted_simulations_data_points[idx] = single_file_col_data / self.initial_number_of_nuclei
                idx += 1

    def calc_p_values_per_frame(self):
        for frame_idx, exp_frame_dist in enumerate(self.exp_data_points_normalized):
            # for random
            frame_data_from_simulations = self.random_simulations_data_points[:, frame_idx]
            self.random_simulations_dist_averages[frame_idx] = np.average(frame_data_from_simulations)
            num_of_equal_or_above = len(frame_data_from_simulations[frame_data_from_simulations >= exp_frame_dist])
            self.p_values_for_random_per_frame[frame_idx] = 0.001 if num_of_equal_or_above / self.number_of_replications <= 0.001 else num_of_equal_or_above / self.number_of_replications
            # for weighted
            frame_data_from_simulations = self.weighted_simulations_data_points[:, frame_idx]
            self.weighted_simulations_dist_averages[frame_idx] = np.average(frame_data_from_simulations)
            num_of_equal_or_above = len(frame_data_from_simulations[frame_data_from_simulations >= exp_frame_dist])
            self.p_values_for_weighted_per_frame[
                frame_idx] = 0.001 if num_of_equal_or_above / self.number_of_replications <= 0.001 else num_of_equal_or_above / self.number_of_replications

    def plot_col_data_with_p_values(self, save_csv: bool = True, from_time_frame: int = 0):
        fig, ax = plt.subplots()
        lines_and_legends = {}
        # plot simulation weighted dist averages
        coding = 'go'
        sim_data = self.weighted_simulations_dist_averages[from_time_frame:]
        t = list(range(from_time_frame, len(self.weighted_simulations_dist_averages)))
        ax.plot(t, sim_data, coding)
        lines_and_legends[
            '{} Simulation Distribution Averages'.format(self.type_of_simulation[1].capitalize())] = coding
        # plot simulation random dist averages
        coding = 'bo'
        sim_data = self.random_simulations_dist_averages[from_time_frame:]
        ax.plot(t, sim_data, coding)
        lines_and_legends['{} Simulation Distribution Averages'.format(self.type_of_simulation[0].capitalize())] = coding
        # plot true exp dist
        coding = 'ro'
        exp_data = self.exp_data_points_normalized[from_time_frame:]
        ax.plot(t, exp_data, coding)
        lines_and_legends['True Experiment Distribution'] = coding
        # plot random p values
        coding = 'k*'
        pval = self.p_values_for_random_per_frame[from_time_frame:]
        ax.plot(t, pval, coding)
        lines_and_legends['RandomVsTrue P-Values<='] = coding
        # plot weighted p values
        coding = 'c*'
        pval = self.p_values_for_weighted_per_frame[from_time_frame:]
        ax.plot(t, pval, coding)
        lines_and_legends['WeightedVsTrue P-Values<='] = coding

        ax.set_title('EXP:{0}\nBootstrap P-Values with each type of policy'.format(self.exp_prefix))
        ax.set_ylabel('fraction')
        ax.set_xlabel('time (h)')

        artist_patches = []
        for leg_label, coding in lines_and_legends.items():
            artist_patches.append(mlines.Line2D([], [],
                                                color=coding[:1],
                                                marker=coding[1:],
                                                label=leg_label))
        # annotate random p_values
        # for idx, p_value in enumerate(self.p_values_for_random_per_frame):
        #     text = '{:.3f}≥'.format(p_value) if p_value != 1 else str(1)
        #     ax.annotate(text, xy=(idx, p_value), xytext=(idx+0.05, p_value), fontsize=6, rotation=40)
        #
        #     # annotate weighted p_values
        #     for idx, p_value in enumerate(self.p_values_for_weighted_per_frame):
        #         text = '{:.3f}≥'.format(p_value) if p_value != 1 else str(1)
        #         ax.annotate(text, xy=(idx, p_value), xytext=(idx + 0.05, p_value), fontsize=6, rotation=100)

        ax.legend(handles=artist_patches, loc=6)
        plt.show()
        path_to_plot = os.sep.join([MAIN_DATA_DIRECTORY_PATH, self.exp_prefix, 'FinalReport',
                                           'Summarized exp: {0}_StartingFrame:{1}.eps'.format(self.exp_prefix, from_time_frame)])
        fig.savefig(path_to_plot, dpi=300)
        #         save to csv
        if save_csv:
            df = pd.DataFrame({
                'Frame Number': ['#{}'.format(x) for x in range(len(self.exp_data_points_normalized[from_time_frame:]))],
                'Fraction of N≥4 in Experiment': self.exp_data_points_normalized[from_time_frame:],
                'Average Fraction of N≥4 in Random Simulations': self.random_simulations_dist_averages[from_time_frame:],
                'Average Fraction of N≥4 in Weighted Simulations': self.weighted_simulations_dist_averages[from_time_frame:],
                'P-Values Random': self.p_values_for_random_per_frame[from_time_frame:],
                'P-Values Weighted': self.p_values_for_weighted_per_frame[from_time_frame:]
            })
            df.to_csv(os.sep.join([MAIN_DATA_DIRECTORY_PATH, self.exp_prefix, 'FinalReport',
                                   'Summarized exp: {0}_StartingFrame:{1}.csv'.format(self.exp_prefix, from_time_frame)]))


for exp_prefix in ['200203_S17', '200203_S19', '200203_S22', '200203_S24', '200604_S06', '200604_S09']:
    # exp_prefix = '200604_S09'
    type_of_simulation = ['randoms', 'weighted']
    col_to_calc_by = 'CTR4total'
    number_of_replications = 1000
    frame_num = 16

    bspe = BootstrappingPerExperiment(exp_prefix, type_of_simulation, col_to_calc_by, number_of_replications, frame_num)
    bspe.read_simulations_col_data_points()
    bspe.calc_p_values_per_frame()
    bspe.plot_col_data_with_p_values(from_time_frame=0)
