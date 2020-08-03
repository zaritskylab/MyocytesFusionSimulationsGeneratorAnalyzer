import os
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.lines as mlines


font = {'size': 8}

matplotlib.rc('font', **font)


def autolabel(rects, ax):
    for rect in rects:
        height = rect.get_height()
        annotation = '{:.2f}'.format(height)

        ax.annotate(annotation,
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')


def plot_exp_vs_simulations(exp_directory_path, ctrs_to_plot: [str]):
    exp_prefix = ' '.join(exp_directory_path.split(os.sep)[-1].split('.')[0].split('_')[:2])
    files = (file for file in os.listdir(exp_directory_path)
             if os.path.isfile(os.path.join(exp_directory_path, file)))
    fig, ax = plt.subplots()
    context_dfs = {}
    raw_colors = ['b', 'g', 'r', 'c']
    colors = {}
    for idx, ctr in enumerate(ctrs_to_plot):
        colors[ctr] = raw_colors[idx]
    marker_types = {
        'true': 'P',
        'simR': 'x',
        'simW': 'd'
    }
    for filename in files:
        if filename != '.DS_Store' and filename.endswith('.csv'):
            exp_path = os.sep.join([exp_directory_path, filename])
            context_dfs[filename.split('.')[0]] = pd.read_csv(exp_path)
    lines_and_legends = {}
    for idx, ctr_to_plot in enumerate(ctrs_to_plot):
        for key, df in context_dfs.items():
            # sum_of_totals = df.loc[:, 'CTRn=1'].values + \
            #                 df.loc[:, 'CTRn=2'].values + \
            #                 df.loc[:, 'CTRn=3'].values + \
            #                 df.loc[:, 'CTRn=4'].values
            initial_number_of_nuclei = np.max(df.loc[:, 'CTRn=1'].values)
            coding = ''
            # coding = colors[ctr_to_plot]
            policy = ''
            if key.find('weighted') != -1:
                coding += 'b'
                coding += marker_types['simW']
                policy = 'weighted by size merge'
            elif key.find('random') != -1:
                coding += 'g'
                coding += marker_types['simR']
                policy = 'random merge'
            else:
                coding += 'r'
                coding += marker_types['true']
                policy = 'true experiment'
            # with 1
            normalized = df.loc[:, ctr_to_plot].values/initial_number_of_nuclei
            # without 1
            if ctr_to_plot == 'CTRn=1':
                continue
            # normalized = df.loc[:, ['CTRn=2', 'CTRn=3', 'CTRn=4']]
            ax.plot(normalized, coding)
            ax.set_xticks(np.arange(len(normalized)))
            ax.set_title(exp_prefix)
            ax.set_ylabel('fraction')
            ax.set_xlabel('time (h)')
            lines_and_legends['{0},{1}'.format(policy, ctr_to_plot)] = coding
            # ax.legend()
    artist_patches = []
    for leg_label, coding in lines_and_legends.items():
        artist_patches.append(mlines.Line2D([], [],
                                            color=coding[:1],
                                            marker=coding[1:],
                                            label=leg_label))
    ax.legend(handles=artist_patches, loc=2)
    plt.show()
    # fig.savefig('/Users/yishaiazabary/Desktop/My Career/OriAvinoamSimulation/Data/{}/{}only_n=4.png'.format('_'.join(exp_prefix.split(' ')),'_'.join(exp_prefix.split(' '))))


def plot_heatmaps(n1_values, n2_values, n3_values, n4_values, enable_log, exp_prefix, frames_labels):
    all_values = [n1_values, n2_values, n3_values, n4_values]
    bytime_and_size_array = np.ndarray((len(all_values), len(frames_labels)))
    for idx, count_bin in enumerate(bytime_and_size_array):
        bytime_and_size_array[idx] = all_values[idx]

    fig, ax = plt.subplots()
    bylog = np.log(bytime_and_size_array)
    # bylog = np.where(bylog==np.NINF, 0, bylog)
    im = ax.imshow(bylog) #, vmin=0, vmax=1)
    ax.set_xticks(np.arange(len(frames_labels)))
    ax.set_yticks(np.arange(len(['n=1', 'n=2', 'n=3', 'n>=4'])))
    ax.set_yticklabels(['n=1', 'n=2', 'n=3', 'n>=4'])
    ax.set_xticklabels(['{}'.format(x) for x in range(len(frames_labels))])
    ax.set_title("{}".format(exp_prefix))
    ax.set_ylabel('Bin')
    ax.set_xlabel('Frame #')
    fig.tight_layout()
    # plt.show()
    fig.savefig('/Users/yishaiazabary/Desktop/My Career/OriAvinoamSimulation/Data/{}/heatmaps/{}.png'.format(exp_prefix,exp_prefix))


def plot_heatmaps_without_n1(n1_values, n2_values, n3_values, n4_values, enable_log, exp_prefix, frames_labels):
    all_values = [n2_values, n3_values, n4_values]
    bytime_and_size_array = np.ndarray((len(all_values), len(frames_labels)))
    for idx, count_bin in enumerate(bytime_and_size_array):
        bytime_and_size_array[idx] = all_values[idx]

    fig, ax = plt.subplots()
    bylog = np.log(bytime_and_size_array)
    # bylog = np.where(bylog==np.NINF, 0, bylog)
    im = ax.imshow(bylog) #, vmin=0, vmax=1)
    ax.set_xticks(np.arange(len(frames_labels)))
    ax.set_yticks(np.arange(len(['n=2', 'n=3', 'n>=4'])))
    ax.set_yticklabels([ 'n=2', 'n=3', 'n>=4'])
    ax.set_xticklabels(['{}'.format(x) for x in range(len(frames_labels))])
    ax.set_title("{}".format(exp_prefix))
    ax.set_ylabel('Bin')
    ax.set_xlabel('Frame #')
    fig.tight_layout()
    plt.show()
    # fig.savefig('/Users/yishaiazabary/Desktop/My Career/OriAvinoamSimulation/Data/{}/heatmaps/{}.png'.format(exp_prefix,exp_prefix))


def plot_counts_bar(n1_values,
                    n2_values,
                    n3_values,
                    n4_values,
                    enable_log,
                    cell_nuclei: str,
                    exp_prefix: str,
                    frames_labels,
                    by_which):
    x = np.arange(len(frames_labels))  # the label locations
    width = 0.15  # the width of the bars

    fig, ax = plt.subplots()
    rects1 = ax.bar(x - ((width / 2) * 2), n1_values, width, label='N=1')
    rects2 = ax.bar(x - width / 2, n2_values, width, label='N=2')
    rects3 = ax.bar(x + width / 2, n3_values, width, label='N=3')
    rects4 = ax.bar(x + ((width / 2) * 2), n4_values, width, label='N>=4')

    # rects1/ = ax.bar(x, n1_values, width, label='N=1')
    # rects2 = ax.bar(x, n2_values, width, label='N=2')
    # rects3 = ax.bar(x, n3_values, width, label='N=3')
    # rects4 = ax.bar(x, n4_values, width, label='N=4')

    if enable_log:
        ax.set_ylabel('Log of Count')
        ax.set_title('{0}: Log(Total {1} fraction) Per Frame'.format(exp_prefix, cell_nuclei))
    else:
        ax.set_ylabel('Count')
        ax.set_title('{0}: {1} fraction Per Frame'.format(exp_prefix, cell_nuclei))

    ax.set_xticks(x)
    ax.set_xticklabels(frames_labels)
    ax.legend()

    autolabel(rects1, ax)
    autolabel(rects2, ax)
    autolabel(rects3, ax)
    autolabel(rects4, ax)
    fig.tight_layout()

    plt.show()
    # fig.savefig('Plots/Data/{0}_{1}.png'.format(exp_prefix, by_which))


def calculations_by_cells_or_nuclei(by_which: str, exp_prefix, exp_df, frame_labels):
    # not normalized
    cntc = []
    if by_which is 'cell':
        cntc = ['CTRn=1', 'CTRn=2', 'CTRn=3', 'CTRn=4']
    else:
        cntc = ['CTR1total', 'CTR2total', 'CTR3total', 'CTR4total']
    # if enable_log:
    #     n1_values = np.log(exp_df.loc[:, cntc[0]].values)
    #     n2_values = np.log(exp_df.loc[:, cntc[1]].values)
    #     n3_values = np.log(exp_df.loc[:, cntc[2]].values)
    #     n4_values = np.log(exp_df.loc[:, cntc[3]].values)
    # else:
    #     n1_values = exp_df.loc[:, cntc[0]].values
    #     n2_values = exp_df.loc[:, cntc[1]].values
    #     n3_values = exp_df.loc[:, cntc[2]].values
    #     n4_values = exp_df.loc[:, cntc[3]].values
    # not normalized
    # plot_counts_bar(n1_values,
    #                 n2_values,
    #                 n3_values,
    #                 n4_values,
    #                 enable_log, cell_nuclei=by_which, exp_prefix)
    # normalized
    sum_of_totals = exp_df.loc[:, cntc[0]].values + \
                    exp_df.loc[:, cntc[1]].values + \
                    exp_df.loc[:, cntc[2]].values + \
                    exp_df.loc[:, cntc[3]].values
    n1_values = exp_df.loc[:, cntc[0]].values / sum_of_totals
    n2_values = exp_df.loc[:, cntc[1]].values / sum_of_totals
    n3_values = exp_df.loc[:, cntc[2]].values / sum_of_totals
    n4_values = exp_df.loc[:, cntc[3]].values / sum_of_totals
    plot_heatmaps_without_n1(n1_values,
                  n2_values,
                  n3_values,
                  n4_values,
                  False,
                  exp_prefix=exp_prefix,
                  frames_labels=frame_labels)
    # plot_heatmaps(n1_values,
    #               n2_values,
    #               n3_values,
    #               n4_values,
    #               False,
    #               exp_prefix=exp_prefix,
    #               frames_labels=frame_labels)
    # plot_counts_bar(n1_values,
    #                 n2_values,
    #                 n3_values,
    #                 n4_values,
    #                 False,
    #                 cell_nuclei=by_which,
    #                 exp_prefix=exp_prefix,
    #                 frames_labels=frame_labels,
    #                 by_which=by_which)


def single_file(
        exp_path: str = '/Users/yishaiazabary/Desktop/My Career/OriAvinoamSimulation/Data/200203_S17_processed.csv',
        frames_labels=[]):
    exp_prefix = exp_path.split(os.sep)[-1].split('.')[0]
    exp_df = pd.read_csv(exp_path)
    frames_labels = ['F.#{}'.format(frame_num) for frame_num in exp_df.loc[:, 'Frame#'].values]

    calculations_by_cells_or_nuclei('cell', exp_prefix, exp_df, frames_labels)
    # calculations_by_cells_or_nuclei('nuclei', exp_prefix, exp_df, frames_labels)


def all_directory(directory_path: str):
    files = (file for file in os.listdir(directory_path)
             if os.path.isfile(os.path.join(directory_path, file)))
    for filename in files:
        if filename != '.DS_Store':
            exp_path = os.sep.join([directory_path, filename])
            single_file(exp_path)


# all_directory('/Users/yishaiazabary/Desktop/My Career/OriAvinoamSimulation/Data/200203_S17')
plot_exp_vs_simulations(exp_directory_path='/Users/yishaiazabary/PycharmProjects/MucsleCellFusionSimulation/Data/200203_S17',
                        ctrs_to_plot=['CTR4total'])