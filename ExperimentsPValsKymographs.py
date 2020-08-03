import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, Colormap


MAIN_DIR = '/Users/yishaiazabary/PycharmProjects/MucsleCellFusionSimulation/Data'


def create_kymograph(list_of_exps_prefixes,
                     starting_frame: int = 0,
                     total_number_of_frames: int = 16):
    all_files_randoms_p_vals = np.zeros((len(list_of_exps_prefixes), total_number_of_frames))
    all_files_weighted_p_vals = np.zeros((len(list_of_exps_prefixes), total_number_of_frames))
    for exp_idx, exp_name in enumerate(list_of_exps_prefixes):
        path_to_exp_summarized = os.sep.join([MAIN_DIR, exp_name, 'FinalReport', 'Summarized exp: {0}_StartingFrame:{1}.csv'.format(exp_name, starting_frame)])
        exp_df = pd.read_csv(path_to_exp_summarized)
        all_files_randoms_p_vals[exp_idx, :] = exp_df.loc[:, 'P-Values Random']
        all_files_weighted_p_vals[exp_idx, :] = exp_df.loc[:, 'P-Values Weighted']

    fig, axis = plt.subplots(2)
    randoms_kymograph = axis[0].imshow(all_files_randoms_p_vals, vmax=0.05, vmin=0, cmap='summer')
    axis[0].set_title('Random Simulations P-values\n for each experiment')
    axis[0].set_xticks(np.arange(total_number_of_frames))
    axis[0].set_yticks(np.arange(len(list_of_exps_prefixes)))
    axis[0].set_xticklabels(['#{}'.format(x) for x in range(total_number_of_frames)], rotation=45)
    # axis[0].set_yticklabels(['Exp No. {}'.format(x) for x in range(len(list_of_exps_prefixes))])
    axis[0].set_yticklabels(list_of_exps_prefixes)
    axis[0].xlabel = 'Time frame (h)'
    axis[0].ylabel = 'Experiment'
    cbar = axis[0].figure.colorbar(randoms_kymograph, ax=axis[0])
    cbar.ax.set_ylabel('P-Value', rotation=-90, va="bottom")
    weighted_kymograph = axis[1].imshow(all_files_weighted_p_vals, vmax=0.05, vmin=0, cmap='summer')

    axis[1].set_title('Weighted Simulations P-values\n for each experiment')
    axis[1].set_xticks(np.arange(total_number_of_frames))
    axis[1].set_yticks(np.arange(len(list_of_exps_prefixes)))
    axis[1].set_xticklabels(['#{}'.format(x) for x in range(total_number_of_frames)], rotation=45)
    # axis[1].set_yticklabels(['Exp No. {}'.format(x) for x in range(len(list_of_exps_prefixes))])
    axis[1].set_yticklabels(list_of_exps_prefixes)
    axis[1].xlabel = 'Time frame (h)'
    axis[1].ylabel = 'Experiment'
    cbar = axis[1].figure.colorbar(randoms_kymograph, ax=axis[1])
    cbar.ax.set_ylabel('P-Value', rotation=-90, va="bottom")

    fig.tight_layout()
    plt.show()
    fig.savefig(os.sep.join([MAIN_DIR, 'pValsKymographs.png']), dpi=300)


example_experiments_list = ['200203_S17', '200203_S19', '200203_S22', '200203_S24', '200604_S06', '200604_S09']
create_kymograph(list_of_exps_prefixes=example_experiments_list)
