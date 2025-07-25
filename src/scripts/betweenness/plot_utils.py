import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import os

from src.utils.plot_utils import *

algo_map = {'genemaniaplus':'RL', 'rwr': 'RWR'}
# net_map = {'STRING11-5-700':'STRING', 'BioGRID-Physical': 'BioGRID-Phy', 'BioGRID-Y2H':'BioGRID-Y2H'}
net_map = {'STRING11-5-700':'STRING', 'BioGRID-Physical': 'BioGRID-Phy', 'BioGRID-Y2H':'BioGRID-Y2H', 'Signor':'Signor'}

def plot_hypergeom_pval(all_criteria_overlap_pvals_topks, interesting_in_pos, interesting_in_top,
                        interesting_in_net, title, filename, ks=[],
                        rank_criteria=['betweenness']):
    '''
        frac_overlap_hypergeom_pvals_topks is a dict. Where, k=each k in topks
        and the value is a tuple(x,y)=> x is fraction of overlapping prots in topk and interesting prot,
        y is the pvalue of that overlap.
    '''

    linestyle_dict = {'betweenness':'solid', 'contr_pathlen_2':'dotted',
                        'contr_pathlen_3':(0,(1,10)),'contr_pathlen_4':'dashdot', 'score':'solid'}
    a_sig = 0.01

    fig, ax = plt.subplots(1)
    for rank_criteron in rank_criteria:
        frac_overlap_hypergeom_pvals_topks = all_criteria_overlap_pvals_topks[rank_criteron]
        # x = list(frac_overlap_hypergeom_pvals_topks.keys())
        # y = [frac_overlap for (frac_overlap, pval) in list(frac_overlap_hypergeom_pvals_topks.values())]
        # pvals = [pval for (frac_overlap, pval) in list(frac_overlap_hypergeom_pvals_topks.values())]

        x = list(frac_overlap_hypergeom_pvals_topks.keys())
        #now keep only the ks present in passed ks
        x = [i for i in x if i in ks]
        y = [frac_overlap for (frac_overlap, pval) in [frac_overlap_hypergeom_pvals_topks[key] for key in x]]
        pvals = [pval for (frac_overlap, pval) in [frac_overlap_hypergeom_pvals_topks[key] for key in x]]
        c = ['g' if i<a_sig else 'r' for i in pvals]

        lines = [((x0, y0), (x1, y1)) for x0, y0, x1, y1 in zip(x[:-1], y[:-1], x[1:], y[1:])]
        colored_lines = LineCollection(lines, colors=c, linewidths=(2,), linestyle=linestyle_dict[rank_criteron],
                                       label=rank_criteron)

        # plot data
        ax.add_collection(colored_lines)
        ax.autoscale_view()

    plt.axhline(y=interesting_in_pos, color='c', linestyle='--', label = 'frac_in_positive')
    plt.axhline(y=interesting_in_top, color='m', linestyle='--', label='frac_in_top_preds')
    plt.axhline(y=interesting_in_net, color='b', linestyle='--', label='frac_in_net')

    legend = ax.legend(loc='upper right', shadow=True, fontsize='x-large',
                       fancybox=True, framealpha=0.5)
    plt.ylim([0, 1])
    plt.xlabel('rank')
    plt.ylabel('fraction of overlap')
    plt.title(title)
    plt.tight_layout()

    # filename1 = filename.replace('.pdf','_'+rank_criteron+'.pdf')
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    plt.savefig(filename)
    plt.savefig(filename.replace('.pdf', '.png'))  # save plot in .png format as well
    plt.close()
    print('Save overlap fig to ', filename)

def plot_hypergeom_pval_multinet(all_criteria_overlap_pvals_topks_multinet,
                        interesting_in_multinet, title=None, filename=None, ks=[]):
    '''
        all_criteria_overlap_pvals_topks_multinet is a dict of dict.
        Where, first key=network_name, second key=k=each k in topks
        and the value is a tuple(x,y)=> x is fraction of overlapping prots in topk and interesting prot,
        y is the pvalue of that overlap.
        Plot overlap between top ranked betweenness proteins and prots of interest.
    '''

    linestyle_dict = {'STRING11-5-700': 'solid', 'BioGRID-Physical': 'dotted',
                      'BioGRID-Y2H': 'dashdot' , 'HI-Union': (0, (1, 10))}

    a_sig = 0.01
    fig, ax = plt.subplots(1)
    for net in all_criteria_overlap_pvals_topks_multinet:
        all_criteria_overlap_pvals_topks = all_criteria_overlap_pvals_topks_multinet[net]
        frac_overlap_hypergeom_pvals_topks = all_criteria_overlap_pvals_topks['betweenness']
        x = list(frac_overlap_hypergeom_pvals_topks.keys())
        #now keep only the ks present in passed ks
        x = [i for i in x if i in ks]
        y = [frac_overlap for (frac_overlap, pval) in [frac_overlap_hypergeom_pvals_topks[key] for key in x]]
        pvals = [pval for (frac_overlap, pval) in [frac_overlap_hypergeom_pvals_topks[key] for key in x]]
        c = ['g' if i<a_sig else 'r' for i in pvals]

        lines = [((x0, y0), (x1, y1)) for x0, y0, x1, y1 in zip(x[:-1], y[:-1], x[1:], y[1:])]
        colored_lines = LineCollection(lines, colors=c, linewidths=(2,), linestyle=linestyle_dict[net],
                                       label=net_name_alias_btns_overlap[net])

        # plot data
        ax.add_collection(colored_lines)
        ax.autoscale_view()

        plt.axhline(y=interesting_in_multinet[net], color='b', linestyle=linestyle_dict[net])

    legend = ax.legend(loc='upper right', shadow=True, fontsize='x-large',
                       fancybox=True, framealpha=0.5)
    plt.ylim([0, 1])
    plt.xlabel('rank')
    plt.ylabel('fraction of overlap')
    plt.title(title)
    plt.tight_layout()

    # filename1 = filename.replace('.pdf','_'+rank_criteron+'.pdf')
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    plt.savefig(filename)
    plt.savefig(filename.replace('.pdf', '.png'))  # save plot in .png format as well
    plt.close()
    print('Save overlap fig to ', filename)



def heatmap_hypergeom_pval_multialg_multinet(all_criteria_overlap_pvals_topks_multialg_multinet,filename=None):

    # Convert to DataFrame
    overlap_df = pd.DataFrame([
        [outer_key, inner_1_key,inner_3_key, value[0],value[1]]
        for outer_key, inner_1_dict in all_criteria_overlap_pvals_topks_multialg_multinet.items()
        for inner_1_key, inner_2_dict in inner_1_dict.items()
        for inner_2_key, inner_3_dict in  inner_2_dict.items()
        for inner_3_key, value in inner_3_dict.items()

    ], columns=['Algorithms', 'Networks', 'top-K-btns', 'Fraction of overlap', 'P-value'])
    overlap_df['Networks'] = overlap_df['Networks'].astype(str).apply(lambda x: net_map.get(x, x))
    algorithms = overlap_df['Algorithms'].unique()
    network_names = overlap_df['Networks'].unique()
    #plot heatmap

    num_algorithms = len(algorithms)
    fig, axes = plt.subplots(1, num_algorithms, figsize=(2.5 * num_algorithms, 5), sharey=True)
    # If there's only one algorithm, ensure 'axes' is treated as a list for consistency
    if num_algorithms == 1:
        axes = [axes]
    # Loop through each algorithm and create a grouped bar plot
    for i, algo in enumerate(algorithms):
        algo_df = overlap_df[overlap_df['Algorithms'] == algo]
        heatmap_data = algo_df.pivot(index='top-K-btns', columns='Networks', values='Fraction of overlap')
        #sort networks
        ordered_netnames = [x for x in net_map.values() if x in network_names ]
        # heatmap_data=heatmap_data[list(net_map.values())]
        heatmap_data=heatmap_data[ordered_netnames]


        orig_min = heatmap_data.min().min()
        orig_max = heatmap_data.max().max()

        if (orig_max-orig_min)<0.3:
            val = 3
            step_size = 0.03
        else:
            val = 5  # wanna keep ticks like 0.24, 0.3 i.e., divisible by 5
            step_size = 0.05

        min_value = int(int(orig_min*100)/val)*val/100
        max_value = int(int(orig_max*100)/val)*val/100+0.00001

        ticks = np.arange(min_value, max_value, step_size)
        sns.heatmap(heatmap_data, cmap="viridis", cbar=True, ax=axes[i],cbar_kws={"ticks": ticks},  linewidths=0.2, linecolor='white')

        axes[i].set_title(alg_plot_name[algo])
        axes[i].set_xlabel('Algorithms', fontsize=12)
        axes[i].set_ylabel('Rank based on betweenness' if i == 0 else '', fontsize=12)  # Only label y-axis for the first plot
        axes[i].set_xticklabels(axes[i].get_xticklabels(), fontsize=10, rotation=45, ha='center')  # Adjust the fontsize as needed
        if i==0:
            axes[i].set_yticklabels(axes[i].get_yticklabels(), fontsize=10)  # Adjust the fontsize as needed

    # Adjust layout
    filename = filename+'heatmap.pdf'
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    plt.tight_layout()
    plt.savefig(filename, bbox_inches='tight')
    plt.savefig(filename.replace('.pdf', '.png'), bbox_inches='tight')  # save plot in .png format as well
    plt.close()
    print('Save overlap fig to ', filename)


def barplot_hypergeom_pval_multialg_multinet(all_criteria_overlap_pvals_topks_multialg_multinet,
                        interesting_in_multialg_multinet,filename=None, ks=[]):

    # Convert to DataFrame
    overlap_df = pd.DataFrame([
        [outer_key, inner_1_key,inner_3_key, value[0],value[1]]
        for outer_key, inner_1_dict in all_criteria_overlap_pvals_topks_multialg_multinet.items()
        for inner_1_key, inner_2_dict in inner_1_dict.items()
        for inner_2_key, inner_3_dict in  inner_2_dict.items()
        for inner_3_key, value in inner_3_dict.items()

    ], columns=['Algorithms', 'Networks', 'top-K-btns', 'Fraction of overlap', 'P-value'])

    #plot barplots
    algorithms = overlap_df['Algorithms'].unique()
    nets = overlap_df['Networks'].unique()
    # Create a color palette for networks
    network_palette = sns.color_palette('Set2', len(nets))
    network_colors = dict(zip(nets, network_palette))

    num_algorithms = len(algorithms)
    fig, axes = plt.subplots(1, num_algorithms, figsize=(6 * num_algorithms, 6), sharey=True)
    # If there's only one algorithm, ensure 'axes' is treated as a list for consistency
    if num_algorithms == 1:
        axes = [axes]
    # Loop through each algorithm and create a grouped bar plot
    for i, algo in enumerate(algorithms):
        algo_df = overlap_df[overlap_df['Algorithms'] == algo]
        sns.barplot(x='top-K-btns', y='Fraction of overlap', hue='Networks', data=algo_df, palette=network_colors,ax=axes[i])
        axes[i].set_title(alg_plot_name[algo])
        axes[i].set_xlabel('Ranks')
        axes[i].set_ylabel('Fraction of overlap' if i == 0 else '')  # Only label y-axis for the first plot


        #plot lines
        for net in nets:
            axes[i].axhline(y=interesting_in_multialg_multinet[algo][net], color=network_colors[net])

    # Adjust layout
    filename = filename+'barplot.pdf'
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    plt.tight_layout()
    plt.legend(loc='upper right', bbox_to_anchor=(1.15, 1))
    plt.savefig(filename, bbox_inches='tight')
    plt.savefig(filename.replace('.pdf', '.png'), bbox_inches='tight')  # save plot in .png format as well
    plt.close()
    print('Save overlap fig to ', filename)




def lineplot_hypergeom_pval_multialg_multinet(all_criteria_overlap_pvals_topks_multialg_multinet,
                                              interesting_in_multialg_multinet, filename=None, ks=[]):
    '''
        all_criteria_overlap_pvals_topks_multialg_multinet is a dict of dict of dict.
        Where, fisrt key=alg_name, second key=network name, third key=k=each k in topks
        and the value is a tuple(x,y)=> x is fraction of overlapping prots in topk and interesting prot,
        y is the pvalue of that overlap.
        Plot overlap between top ranked betweenness proteins and prots of interest.
    '''

    linestyle_dict = {'STRING11-5-700': 'solid', 'BioGRID-Physical': 'dotted',
                      'BioGRID-Y2H': 'dashdot' , 'HI-Union': (0, (1, 10)), 'Signor': 'solid'}

    a_sig = 0.01
    n_algs = len(all_criteria_overlap_pvals_topks_multialg_multinet.keys())
    if n_algs==1:
        fig, ax = plt.subplots(1)
        axs=[ax]

    elif n_algs==2:
        fig, (ax1, ax2) = plt.subplots(1,2)
        axs=[ax1,ax2]
    count=0
    for alg in all_criteria_overlap_pvals_topks_multialg_multinet:
        all_criteria_overlap_pvals_topks_multinet = all_criteria_overlap_pvals_topks_multialg_multinet[alg]
        interesting_in_multinet = interesting_in_multialg_multinet[alg]
        for net in all_criteria_overlap_pvals_topks_multinet:
            all_criteria_overlap_pvals_topks = all_criteria_overlap_pvals_topks_multinet[net]
            frac_overlap_hypergeom_pvals_topks = all_criteria_overlap_pvals_topks['betweenness']
            x = list(frac_overlap_hypergeom_pvals_topks.keys())
            #now keep only the ks present in passed ks
            x = [i for i in x if i in ks]
            y = [frac_overlap for (frac_overlap, pval) in [frac_overlap_hypergeom_pvals_topks[key] for key in x]]
            pvals = [pval for (frac_overlap, pval) in [frac_overlap_hypergeom_pvals_topks[key] for key in x]]
            c = ['g' if i<a_sig else 'r' for i in pvals]

            lines = [((x0, y0), (x1, y1)) for x0, y0, x1, y1 in zip(x[:-1], y[:-1], x[1:], y[1:])]
            colored_lines = LineCollection(lines, colors=c, linewidths=(2,), linestyle=linestyle_dict[net],
                                           label=net_name_alias_btns_overlap.get(net,net))

            # plot data
            axs[count].add_collection(colored_lines)
            axs[count].autoscale_view()

            axs[count].axhline(y=interesting_in_multinet[net], color='b', linestyle=linestyle_dict[net])
            axs[count].set_ylim(0,1)
            axs[count].set_xlabel('ranks')
            axs[count].set_ylabel('fraction of overlap')
            axs[count].set_title(alg_plot_name[alg])
            axs[count].set_box_aspect(0.75)
            legend = axs[count].legend(loc='upper right', shadow=True, fancybox=True)
        count += 1



    # filename1 = filename.replace('.pdf','_'+rank_criteron+'.pdf')
    plt.tight_layout()
    filename = filename+'lineplot.pdf'
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    plt.savefig(filename, bbox_inches='tight')
    plt.savefig(filename.replace('.pdf', '.png'), bbox_inches='tight')  # save plot in .png format as well
    plt.close()
    print('Save overlap fig to ', filename)

def plot_KS(frac_prots_ge_btns_marker,marker, title, filename):
    markers_dict = {'ppi_prots':'o', 'ess_cell': 'v', 'ess_org': 's', 'viral_sars2--':'P'}

    btns_markers = list(frac_prots_ge_btns_marker['ppi_prots'].keys())
    fig,ax = plt.subplots()
    for prot_type in frac_prots_ge_btns_marker:
        plt.plot(btns_markers, list(frac_prots_ge_btns_marker[prot_type].values()),
                   marker=markers_dict[prot_type], label=prot_type)
        # plt.loglog(btns_markers, list(frac_prots_ge_btns_marker[prot_type].values()),
        #         marker =markers_dict[prot_type] ,label = prot_type)
    #commenting out the log scale
    # ax.set_xscale("log", base=10)
    # ax.set_yscale("log", base=10)

    legend = ax.legend(loc='upper left', shadow=True, fontsize='x-large',
                       fancybox=True, framealpha=0.5)
    # Put a nicer background color on the legend.

    plt.xlabel(marker +' along betweenness score')
    plt.ylabel('fraction of prots')
    plt.title(title)
    plt.tight_layout()
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    plt.savefig(filename)
    plt.savefig(filename.replace('.pdf', '.png'))  # save plot in .png format as well
    plt.close()
    print('Save fig to ', filename)


def plot_prots_appearing_at_each_pathlens(new_prots_each_pathlens, filename):
    #For balancing alpha only
    '''Input: a dict with keys: ['network', 'term', 'alg', 'alpha','new_appearing_prots'], each value is a list.
    Especially: value for the key 'new_appearing_prots' is a list of dicts.
    Where inner dict keys:['pathlen_2', 'pathlen_3','pathlen_4']

    Output: a plot where along x-axis we will have networks and along y axis we will have stacked bar chart.
    In the bar chart, we will have the #of_prots appearing as we go along each path_len.
    '''

    #converting prots_appearing_at_each_pathlens into dataframe
    new_prots_each_pathlens_df = pd.DataFrame(new_prots_each_pathlens)
    new_prots_each_pathlens_df = new_prots_each_pathlens_df[['term','alg','network','pathlen_2', 'pathlen_3','pathlen_4']]
    #for a certain term-alg-alpha combo plot all the networks and  #new_appearing_prots in one plot.
    #the following will return a list of unique (term, alg) tuples present in prots_appearing_at_each_pathlens
    term_alg_alpha_list = list(new_prots_each_pathlens_df.groupby(by=['term', 'alg']).groups.keys())

    for (term, alg) in term_alg_alpha_list:
        df = new_prots_each_pathlens_df[(new_prots_each_pathlens_df['term']==term)&
                                        (new_prots_each_pathlens_df['alg']==alg)]

        df.set_index('network', inplace=True)
        df=df[['pathlen_2', 'pathlen_3','pathlen_4']]
        ax = df.plot.bar(stacked=True)

        legend = ax.legend(loc='upper right', shadow=True, fontsize='x-large',
                           fancybox=True, framealpha=0.5)
        # Put a nicer background color on the legend.

        plt.xlabel('Networks')
        plt.ylabel('New prots appearing at each path length')
        plt.title(term + '_'+ alg )
        plt.tight_layout()

        filename1 =filename.replace('.pdf', '_'+term + '_'+ alg +'.pdf')
        os.makedirs(os.path.dirname(filename1), exist_ok=True)
        plt.savefig(filename1)
        plt.savefig(filename1.replace('.pdf', '.png'))  # save plot in .png format as well
        plt.show()
        plt.close()
        print('Save fig to ', filename1)




def barplot_from_dict(n_essential_prots_per_topk, x_label, y_label,ymax, filename,title=''):
    '''
    n_essential_prots_per_topk: dict where key = k, value: number of essential prot in top k
    '''
    X = list(n_essential_prots_per_topk.keys())
    Y =  list(n_essential_prots_per_topk.values())
    sns.barplot(x=X, y=Y)

    plt.ylim([0,ymax])
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(title)
    plt.tight_layout()

    os.makedirs(os.path.dirname(filename), exist_ok=True)
    plt.savefig(filename)
    plt.savefig(filename.replace('.pdf', '.png'))  # save plot in .png format as well
    # plt.show()
    plt.close()
    print('Save fig to ', filename)

def barplot_from_df(n_interesting_prots_per_topk_df, filename, x, y,ymax, title=''):
    '''
    n_essential_prots_per_topk: dict where key = k, value: number of essential prot in top k
    '''
    sns.barplot(n_interesting_prots_per_topk_df, x = x, y = y, hue='alpha')

    plt.ylim([0, ymax])
    plt.title(title)
    plt.tight_layout()

    os.makedirs(os.path.dirname(filename), exist_ok=True)
    plt.savefig(filename)
    plt.savefig(filename.replace('.pdf', '.png'))  # save plot in .png format as well
    # plt.show()
    plt.close()
    print('Save fig to ', filename)

def scatter_plot(X, Y,  x_label, y_label, ymin=0, ymax=100, title='', filename=''):


    #in this case the last to points are for positive/src_bin and to_predicted_proteins_bin. And we want to show them in
    #different color
    color_non_top_pos = ['blue']*(len(X)-2)
    marker_non_top_pos = '.'
    plt.scatter(X[0:-2],Y[0:-2],color = color_non_top_pos, marker=marker_non_top_pos)

    color_pos = ['red']
    marker_pos='^'
    plt.scatter(X[-2],Y[-2],color = color_pos, marker=marker_pos)


    color_top = ['yellow']
    marker_top='s'
    plt.scatter(X[-1],Y[-1],color = color_top, marker=marker_top)



    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.ylim([ymin, ymax])
    plt.title(title)
    plt.tight_layout()

    os.makedirs(os.path.dirname(filename), exist_ok=True)
    plt.savefig(filename)
    plt.savefig(filename.replace('.pdf', '.png'))  # save plot in .png format as well
    # plt.show()
    plt.close()
    print('Save fig to ', filename)

def box_plot(data, x, y, ymin, ymax, title, filename):
    sns.boxplot(data, x=x, y=y )

    plt.ylim([ymin, ymax])
    plt.title(title)
    plt.tight_layout()

    os.makedirs(os.path.dirname(filename), exist_ok=True)
    plt.savefig(filename)
    plt.savefig(filename.replace('.pdf', '.png'))  # save plot in .png format as well
    # plt.show()
    plt.close()
    print('Save fig to ', filename)

def plot_non_zero_contr_nodes(non_zero_contr_df, output_dir):
    plt.figure(figsize=(4, 4))

    non_zero_contr_df['alg_name'] = non_zero_contr_df['alg_name'].astype(str).apply(lambda x: algo_map[x])
    non_zero_contr_df['dataset_name'] = non_zero_contr_df['dataset_name'].astype(str).apply(lambda x: net_map.get(x,x))

    sns.barplot(data=non_zero_contr_df, x='alg_name', y='frac_non_zero_contr_nodes', hue='dataset_name')
    plt.xlabel("Algorithm", fontsize=12)
    plt.ylabel("Fraction of non-zero contributing nodes", fontsize=12)

    # Show the plot
    plt.tight_layout()
    plt.savefig(f'{output_dir}/summary_non_zero_contr_nodes_barplot.pdf')
    plt.show()
    plt.clf()
