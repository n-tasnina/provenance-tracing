import numpy as np
import scipy.stats as stats
import os
import pandas as pd

from src.FastSinkSource.src.algorithms import rl_genemania as gm
from src.FastSinkSource.src.algorithms import PageRank_Matrix as rwr

import rbo
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import gridspec


alg_alias = {'rwr': rwr, 'genemaniaplus': gm}

#************************************* NETWORK PROCESSING *********************************************
def log_inv_mat(M_pathmtx):
    '''
    This will convert each value in the input matrix to 10^-(value) which is equvalent to taking inverse
    of 10 base log.
    Also, it will convert any value=1 in the output matrix to be zero.
    '''
    M = np.power(np.full_like(M_pathmtx, 10), (-1) * M_pathmtx)
    # in M_pathmtx an index with value 0 means no edge. The above statement turns all 0 values to 1.
    # Now to preserve the same meaning for such indices, we need to convert all 1 to 0 in M.
    where_1 = np.where(M == 1)  # TODO: replace 1's with 0's in time efficient way
    M[where_1] = 0
    return M


def get_M_pathmtx_loginv(net_obj, alg_name, alpha):
    fluid_flow_mat_file_M = "%s/fluid-flow-mat-%s-a%s.npy" % (net_obj.out_pref, alg_name,
                                                              str(alpha).replace('.', '_'))
    ##########CREAE or LOAD THE DIFFUSION MAYTRICES
    force_matrix = False
    if (os.path.exists(fluid_flow_mat_file_M) and force_matrix == False):
        M_pathmtx = np.load(fluid_flow_mat_file_M)
        M_pathmtx_loginv = log_inv_mat(M_pathmtx)
    else:
        M_pathmtx_loginv = alg_alias[alg_name].get_M(net_obj.W, alpha)

    return M_pathmtx_loginv

#***************************************** ALGO PREDICTIONS SCORE *********************************************
def get_top_k_predictions(pred_file, alpha, k, orig_pos):
    '''
    This will return a dataframe with proteins sorted descendingly according to their predicted score.
    It removes the original positive proteins from the output dataframe.
    '''

    if not os.path.isfile(pred_file):
        print("Warning: %s not found. skipping" % (pred_file))
    else:
        print("reading %s for alpha=%s" % (pred_file, alpha))
        df = pd.read_csv(pred_file, sep='\t')

        # remove the original positives for downstream analysis
        df = df[~(df['prot'].isin(orig_pos))]
        df.reset_index(inplace=True, drop=True)
        df = df[:k]
        return df
    return None

def term_based_pred_file(pred_file, term):
    pred_file = \
        pred_file.split('.')[0] + '-' +term + '.' + pred_file.split('.')[-1]
    #also in GO term we might see ':' which we replaced by '-' while writing the pediction to file. so consider that here.
    pred_file  = pred_file.replace(':','-')
    return pred_file


def get_balancing_alpha(config_map, dataset, alg_name,  term):
    alpha_summary_filename = config_map['output_settings']['output_dir'] + \
                             "/viz/%s/%s/param_select/" % (dataset['net_version'], dataset[
        'exp_name']) + '/' + alg_name + '/alpha_summary.tsv'
    alpha_summary_df = pd.read_csv(alpha_summary_filename, sep='\t', index_col=None)[
        ['term', 'balancing_alpha']]
    term_2_balancing_alpha_dict = dict(
        zip(alpha_summary_df['term'], alpha_summary_df['balancing_alpha']))

    balancing_alpha = term_2_balancing_alpha_dict[term]
    return balancing_alpha

#************************* NETWORK ANALYSIS*************************
def is_neighbor(W, source, target):
    ''' This will return True if there is an edge from the source to target'''
    ''' Consider the Weight matrix in a format that along rows we have targets, along columns we have sources'''
    if W[target][source]!=0:
        return True
    return False

def find_in_neighbors(W, target):
    ''' Find the nodes from which there are incoming edges to the target'''
    neighbors = np.where(W[target]>0)[0]
    return list(neighbors)

def is_diag_zero(W):
    if W.diagonal().sum()==0:
        return True
    else:
        return False


#********************** STATISTICAL TEST ******************************
def pearson_correlations(list1, list2):
    #write code for finding pearson correlations coefficient here
    return stats.pearsonr(list1, list2)

def kendal_tau(list1, list2):
    kt = stats.kendalltau(list1, list2)
    return kt.correlation, kt.pvalue

#*********************** LIST COMPARE ***************************************
def compute_jaccard(set1, set2):
    '''
    pass tow set or lists.
    '''
    return len(set(set1).intersection(set(set2)))/ len(set(set1).union(set(set2)))


#*************************** FILE WRITE *********************************
def save_dict(dict1, filename):
    #write key and values of a dict in tab separated way
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    f = open(filename, 'w')
    for key in dict1:
        val = "{:.2e}".format(dict1[key])
        f.write(str(key) + '\t' + val +'\n')
    f.close()


#************** Miscellaneous *********************
def track_rank(path_based_contr, scores_df):
    '''
    path_based_contr: dataframe containing contribution along each length of paths for all target nodes
    pred_scores: contains predicted score for all target nodes
    This function finds if we sum up the contribution along paths of length till 4, if the rank of the nodes remain same as if we consider paths of all lengths.
    '''


    path_based_contr.set_index('target', inplace=True)
    path_based_contr['frac_contr_1_4'] = path_based_contr['frac_pathlen_1']+path_based_contr['frac_pathlen_2']+path_based_contr['frac_pathlen_3']+path_based_contr['frac_pathlen_4']
    scores_df.set_index('prot_idx', inplace=True)
    scores_df = scores_df.join(path_based_contr['frac_contr_1_4'], how='left')

    original_order = scores_df.index.values.tolist()

    #compute order if we considered only till paths of length 4
    path_1_4_rank = (scores_df['score']*scores_df['frac_contr_1_4']).sort_values(ascending=False).index.tolist()

    rbo_overlap_all = round( rbo.RankingSimilarity(original_order, path_1_4_rank).rbo(), 4)
    rbo_overlap_1000 = round( rbo.RankingSimilarity(original_order[0:1000], path_1_4_rank[0:1000]).rbo(), 4)

    print(f'RBO overlap:\n all: {rbo_overlap_all}  top 1000: {rbo_overlap_1000}')
    return rbo_overlap_all, rbo_overlap_1000


def frac_contr_plot(target_wise_frac_contr, target_wise_n_contr, n_targets, out_dir, type):
    '''
    target_wise_frac_contr (dict): each key is a target index. Value=list containing the fraction of score contributed by
    non-zero contributing nodes.
    '''
    if n_targets > 100:
        n = 4  # show every 4 ticks to avoid over crowding
    else:
        n = 1
    plt.clf()
    # bin the fraction of contribution
    bins = [0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1]
    bin_labels = ['0-1', '1-2', '2-3', '3-4', '4-5', '5-10', '10-20', '20-30',
                  '30-40', '40-50', '50-100']
    df = pd.DataFrame()
    for prot in target_wise_frac_contr:
        frac_contr = target_wise_frac_contr[prot]
        binned = pd.cut(frac_contr, bins)

        # Count the frequency of values in each bin
        bin_counts = binned.value_counts()
        # Calculate the fraction of values in each bin, i.e., what fraction of node belong to each bin
        bin_fractions = bin_counts / len(frac_contr)
        target_df = pd.DataFrame(
            {'Target': [prot] * len(bin_labels), 'Contribution range': bin_labels, 'frac_nodes': bin_fractions.values})

        df = pd.concat([df, target_df], axis=0)

    # keep n_targets per heatmap
    n_plots = int(len(list(target_wise_frac_contr.keys())) / n_targets)

    for i in range(n_plots):
        # Create a figure with 2 subplots: one for the heatmap and one for the barplot
        df_small = df[(n_targets * i) * len(bin_labels): (n_targets * i + n_targets) * len(bin_labels)]

        gs = gridspec.GridSpec(2, 1, height_ratios=[0.75, 1])

        # fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 10), gridspec_kw={'height_ratios': [1, 1]})
        fig = plt.figure(figsize=(12, 10))
        ax1 = fig.add_subplot(gs[0])  # First row

        # BAR PLOT
        targets = list(df_small['Target'].unique())
        n_nodes = [target_wise_n_contr[x] for x in targets]
        sns.barplot(x=targets, y=n_nodes, ax=ax1)
        ax1.set_ylabel("# of contributing intermediate nodes", fontsize=14)

        # ax1.set_xticklabels(ax1.get_xticklabels(), rotation=90)
        # Get current x-tick positions and labels
        xticks = ax1.get_xticks()
        xticklabels = ax1.get_xticklabels()

        # Show every 5th x-tick
        n = n
        ax1.set_xticks(xticks[::n])  # Set ticks at every 5th position
        ax1.set_xticklabels(
            [label.get_text() for i, label in enumerate(xticklabels) if i % n == 0],
            rotation=90)  # Set labels for those ticks

        # HEATMAP PLOT
        data = df_small.pivot(columns='Target', index='Contribution range', values='frac_nodes')
        # retain the original order of target and contribution range after pivot.
        data = data.reindex(bin_labels)
        data = data.reindex(columns=list(df_small['Target'].unique()))

        ax2 = fig.add_subplot(gs[1])  # Second row, first column
        # sns.heatmap(data, cmap="icefire", cbar=False, linewidths=0.05, square=True, ax=ax2)

        # Create a new axis for the color bar, and control its position manually
        # for 50 targets at a time
        # cbar_ax = fig.add_axes([0.25, 0.25, 0.5, 0.01])  # [left, bottom, width, height]
        # sns.heatmap(data, cmap="icefire", annot=False, linewidths=0.05, square=True, ax=ax2, cbar=True, cbar_ax=cbar_ax, cbar_kws={"orientation": "horizontal"})

        # for all targets
        # cbar_ax = fig.add_axes([0.25, 0.001, 0.5, 0.01])  # [left, bottom, width, height]
        cbar_ax = fig.add_axes([0.15, 0.05, 0.7, 0.02])
        sns.heatmap(data, cmap="icefire", annot=False, linewidths=0.05, ax=ax2, cbar=True, cbar_ax=cbar_ax,
                    cbar_kws={"orientation": "horizontal", "shrink": 0.5})

        # Add labels and title
        ax2.set_xlabel("Targets")
        ax2.set_ylabel("Contribution ranges (%)")

        # # Adjust layout to remove excess white space
        plt.subplots_adjust(hspace=0.4, top=0.5)

        # Show the plot
        plt.tight_layout()
        plt.savefig(f'{out_dir}_{type}_heatmap_barplot_{i}.pdf')
        plt.show()
        plt.clf()


def plot_target_wise_non_zero_contr_nodes(target_wise_n_contr, out_dir, type):
    '''
    target_wise_frac_contr (dict): each key is a target index. Value=list containing the fraction of score contributed by
    non-zero contributing nodes.
    '''
    targets = list(target_wise_n_contr.keys())
    n_nodes = [target_wise_n_contr[x] for x in targets]
    sns.barplot(x=targets, y=n_nodes)
    plt.ylabel("# of contributing intermediate nodes", fontsize=14)

    # Show the plot
    plt.tight_layout()
    plt.savefig(f'{out_dir}_{type}_barplot.pdf')
    plt.show()
    plt.clf()



def global_net_usage(target_wise_btns_file, scores_df, k, orig_pos):
    '''
    target_wise_btns_file = a pickle file that contains the contribution going via each node (along paths of length 2-4) to each target node. Source=coulmns, target=rows.
    '''
    # Load the object from the pickle file


    out_dir = os.path.dirname(target_wise_btns_file)+'/'
    import pickle
    with open(target_wise_btns_file, 'rb') as f:
        target_wise_btns_mat = pickle.load(f)


    #find top k nodes excluding the original positive prots
    scores_df_non_src = scores_df[~scores_df['prot'].isin(orig_pos)]
    non_src_prot_idx = scores_df_non_src['prot_idx'].tolist()
    top_k = scores_df_non_src[:k]['prot_idx'].tolist()
    idx_to_prot = dict(zip(scores_df['prot_idx'], scores_df['prot']))

    #get the scores of nodes in the order of their index
    score_index_ordered_df = scores_df.sort_values('prot_idx')
    score_index_ordered_df.set_index('prot_idx', inplace=True)
    scores_index_ordered =  np.array(score_index_ordered_df['score'])

    assert score_index_ordered_df.index.tolist() == list(range(target_wise_btns_mat.shape[0])), print('not all index present')
    assert target_wise_btns_mat.shape[0] == scores_index_ordered.shape[0], print('mismatch')


    # find out for each top predited node, how many intermediate nodes have non-zero contribution
    all_non_zero_cont_nodes = set()
    target_wise_n_contr =dict()
    target_wise_n_non_src_contr =dict()


    target_wise_frac_contr = dict()
    target_wise_non_src_frac_contr = dict()


    for idx in top_k:
        prot = idx_to_prot[idx]
        score = scores_index_ordered[idx]
        contr_from_intermediate_nodes = target_wise_btns_mat[idx]
        non_zero_contr_nodes = np.nonzero(contr_from_intermediate_nodes)[0]
        n_non_zero_contr_nodes =len(non_zero_contr_nodes)
        target_wise_n_contr[prot]=n_non_zero_contr_nodes

        non_src_non_zero_contr_nodes = list(set(non_zero_contr_nodes).intersection(non_src_prot_idx))
        n_non_src_non_zero_contr_nodes = len(non_src_non_zero_contr_nodes)
        target_wise_n_non_src_contr[prot] = n_non_src_non_zero_contr_nodes
        # print(f'for top pred: {idx}    #of non-zero contributing nodes: {n_non_zero_contr_nodes}')
        # print(f'#of non-source non-zero contributing nodes: {n_non_src_non_zero_contr_nodes}')

        #compute what fraction of the score is contrbuted via each node
        target_wise_frac_contr[prot] = contr_from_intermediate_nodes[non_zero_contr_nodes]/score

        target_wise_non_src_frac_contr[prot] = contr_from_intermediate_nodes[non_src_non_zero_contr_nodes]/score

        all_non_zero_cont_nodes = all_non_zero_cont_nodes.union(set(non_zero_contr_nodes))
    #plot a heatmap for plotting the fraction of score contributed by the nodes
    n_targets = k
    frac_contr_plot(target_wise_frac_contr, target_wise_n_contr, n_targets, out_dir, type='all')
    frac_contr_plot(target_wise_non_src_frac_contr, target_wise_n_non_src_contr, n_targets, out_dir, type='non_src')

    plot_target_wise_non_zero_contr_nodes(target_wise_n_contr, out_dir, type='all')
    plot_target_wise_non_zero_contr_nodes(target_wise_n_contr, out_dir, type='non_src')



    print('non zero contributing nodes across all top preds: ', len(all_non_zero_cont_nodes) )
    return len(all_non_zero_cont_nodes)



































