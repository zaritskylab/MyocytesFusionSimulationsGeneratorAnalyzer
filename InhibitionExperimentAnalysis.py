import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


DATA_PATH = 'Data/Copy of Fusion_index_Analysis__of_SCH+CaMKII_DONE.csv'
SAVEFIG = True
SHOWFIG = True

org_df = pd.read_csv(DATA_PATH)#, index_col='NucNumber')
print(org_df.head())

# Collecting experiments' data by number of nuclei
num_of_experiments = 6
cols_of_experiments = [3,5,7,9,11,13]
exps_data = {}
for exp_idx in cols_of_experiments:
    exp_series = org_df.iloc[:, exp_idx]
    nuc_1 = exp_series.values[0::4]
    nuc_2 = exp_series.values[1::4]
    nuc_3 = exp_series.values[2::4]
    nuc_4 = exp_series.values[3::4]
    exps_data[exp_idx] = {
        'nuc_1': nuc_1,
        'nuc_2': nuc_2,
        'nuc_3': nuc_3,
        'nuc_4': nuc_4
    }
# plot visuals related input
coding = {
    'nuc_1': 'blacko',
    'nuc_2': 'gv',
    'nuc_3': 'r^',
    'nuc_4': 'tab:blueD'
}
legend_labels = {
    'nuc_1': 'mono-nucleated',
    'nuc_2': 'bi-nucleated',
    'nuc_3': 'tri-nucleated',
    'nuc_4': '≥4 nuclei'
}


fig, ax = plt.subplots()

time_frames = np.arange(8, 24, 1)
# get the maximum count of single nuclei within each experiment for normalization
exps_max = [max(exp_data['nuc_1']) for exp_data in exps_data.values()]

for nuc_idx in range(1, 5):
    # collect all experiments' data regarding cells with n nuclei (according to nuc_idx)
    all_nuc_data = np.asarray([exp_data[f'nuc_{nuc_idx}'] for exp_data in exps_data.values()], dtype='float64')
    # normalize to calculate fraction of nuclei
    for single_exp_idx, single_exp_data in enumerate(all_nuc_data):
        all_nuc_data[single_exp_idx] = (single_exp_data*(nuc_idx))/exps_max[single_exp_idx]
    # calculate the mean and standard deviation
    nuc_avg = all_nuc_data.mean(axis=0)
    nuc_std = all_nuc_data.std(axis=0)

    # plot with appropriate style and labeling
    color = coding[f'nuc_{nuc_idx}'][:-1]
    marker = coding[f'nuc_{nuc_idx}'][-1:]
    label = legend_labels[f'nuc_{nuc_idx}']
    ax.errorbar(time_frames, nuc_avg, label=label, yerr=nuc_std, linewidth=0.2,elinewidth=0.5, fmt='-o', color=color, marker=marker, capsize=5)

# plot style save/show
ax.set_xlabel('Time (hours)', size=15, fontweight='semibold' )
ax.set_ylabel('% of total nuclei', size=15, fontweight='semibold' )

ax.set_xticks(np.arange(8, 25, 2))
ax.set_xticklabels([f'{tick}' for tick in np.arange(8, 25, 2)], size=15, fontweight='semibold')

ax.set_yticks(np.arange(0, 1.1, 0.2))
ax.set_yticklabels([f'{int(tick*100)}' for tick in np.arange(0, 1.1, 0.2)], size=15, fontweight='semibold')

plt.legend(frameon=False, loc='best', prop={'weight':'semibold'})
plt.tight_layout()

if SAVEFIG:
    plt.savefig('test.eps', dpi=300, transparent=True, bbox_inches='tight')
if SHOWFIG:
    plt.show()
