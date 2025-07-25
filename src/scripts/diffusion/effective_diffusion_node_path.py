import os, sys
import yaml
import argparse
import numpy as np
import matplotlib
import networkx as nx
import copy
import time

import subprocess
import logging

logging.basicConfig(filename='diffisuion_eppsteins.log', filemode='a', level=logging.INFO, \
                    format='%(message)s')

HIGH_WEIGHT = 1e6

if __name__ == "__main__":
    # Use this to save files remotely. 
    matplotlib.use('Agg')
import pandas as pd
# sys.path.insert(1, "/data/tasnina/Provenance-Tracing/SARS-CoV-2-network-analysis/")
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../.."))
sys.path.insert(1, project_root)

from src.FastSinkSource.src.main import setup_dataset
from src.FastSinkSource.src.utils import config_utils
from src.FastSinkSource.src.algorithms import alg_utils
from src.FastSinkSource.src.algorithms import rl_genemania as gm
from src.FastSinkSource.src.algorithms import PageRank_Matrix as rwr
import src.scripts.utils as script_utils
alg_alias = {'rwr': rwr, 'genemaniaplus': gm}

def parse_args():
    parser = setup_opts()
    args = parser.parse_args()
    kwargs = vars(args)
    with open(args.config, 'r') as conf:
        config_map = yaml.load(conf, Loader=yaml.FullLoader)
        # config_map = yaml.load(conf)
    # TODO check to make sure the inputs are correct in config_map

    return config_map, kwargs


def setup_opts():
    ## Parse command line args.
    parser = argparse.ArgumentParser(description="Script to analyze the top contributors of each prediction's "
                                                 "diffusion score, as well as the effective diffusion (i.e.,"
                                                 " fraction of diffusion received from non-neighbors)")
    # general parameters
    group = parser.add_argument_group('Main Options')
    group.add_argument('--config', type=str, default="/data/tasnina/Provenance-Tracing/SARS-CoV-2-network-analysis/"
                     "fss_inputs/config_files/provenance/string700_biogrid_physical_biogrid_y2h_s1.yaml"
                       , help="Configuration file used when running FSS. ")

    group.add_argument('--run-algs', type=str, action='append', default=[])

    group.add_argument('--k-to-test', '-k', type=int, action='append', default=[332],
                       help="k-value(s) for which to get the top-k predictions to test. " +
                            "If not specified, will check the config file.")
    group.add_argument('--pos-k', action='store_true', default=True,
                       help="if true get the top-k predictions to test is equal to the number of positive annotations")

    group.add_argument('--n-sp', '-n', type=int, default=1000,
                       help="n-sp is the number of shortest paths to be considered" +
                            "If not specified, will check the config file. Default=20")

    group.add_argument('--max-len', type=int, default=3,
                       help="for each top prediction, for how many top contributing sources we wanna analyse the path" +
                            "Default=20")
    group.add_argument('--stat-sig-cutoff', type=float,
                       help="Cutoff on the node p-value for a node to be considered in the topk. " + \
                            "The p-values should already have been computed with run_eval_algs.py")

    group.add_argument('--force-matrix', action='store_true', default=False,
                       help="Compute the diffusion matrix again")
    group.add_argument('--force-ksp', action='store_true', default=False,
                       help="Force re-running the path diffusion analysis")
    group.add_argument('--force-contr', action='store_true', default=False,
                       help="Force re-running the path diffusion analysis")

    # group.add_argument('--plot-only', action='store_true', default=False,
    #                    help="Force re-running the path diffusion analysis")

    group.add_argument('--balancing-alpha-only', action='store_true', default=True,
                       help="Ignore alpha from config file rather take the alpha value\
                            that balanced the two loss terms in quad loss function for the corresponding\
                            network-term-alg")
    return parser

def add_supersource(G, ss_index, pos_idices):
    edges_ss_to_s = list(zip([ss_index]*len(pos_idices), pos_idices))
    edges_s_to_ss = list(zip(pos_idices,[ss_index]*len(pos_idices)))
    G.add_edges_from(edges_ss_to_s, weight = HIGH_WEIGHT )
    G.add_edges_from(edges_s_to_ss, weight = HIGH_WEIGHT )

def process_java_computed_eppstein_output_for_supersource(path_df):
    # process tha output from java code of Eppsitein's
    path_df['cost'] = path_df['cost'].astype(str).apply(lambda x: x.split(':')[0])
    path_df['cost'] = path_df['cost'].apply(lambda x: float(x)-HIGH_WEIGHT)
    # find the actual cost i.e. cost without log implied.
    # multiplying with (10^HIGH_WEIGHT) at the end to remove the edge cost incorporated by edges from supersource
    path_df['actual_cost'] = path_df['cost'].apply(lambda x: 10 ** (-x))

    # path has been saved as source-intermediate1-intermediate2-target in this form. So counting '-' will
    # give the path length
    #Subtract 1 to remove the edge from supersource
    #TODO Make sure that in path supersource_idx comes up only once.
    path_df['length'] = path_df['path'].astype(str).apply(lambda x: x.count('-')-1)
    return path_df

def compute_path_len_wise_contr_for_supersource(target, path_df, R, pos_nodes_idx, path_length_wise_rate_contr):
    # total_rate_contr_from_s_2_t = R[target][source]  # in R, columns are the sources and rows are the targets
    # paths is a tuple as (path_length, actual_cost)
    paths = tuple(zip(path_df['length'], path_df['actual_cost']))
    for path in paths:
        path_length = path[0]
        path_cost = path[1]
        if path_length in path_length_wise_rate_contr[(target)]:
            path_length_wise_rate_contr[(target)][path_length] += path_cost

    # compute fraction of total score coming via each path length
    for path_length in path_length_wise_rate_contr[(target)]:
        path_length_wise_rate_contr[(target)][path_length] /= sum(R[target, pos_nodes_idx])


def write_path_len_wise_contr_for_supersource(target, score, path_length_wise_rate_contr, filename):
    out_f = open(filename, 'a')
    rounded_score = round(score, 6)
    out_str = str(target) + '\t' + str(rounded_score)
    for path_length in path_length_wise_rate_contr[(target)]:
        rounded_val = round(path_length_wise_rate_contr[(target)][path_length], 6)
        out_str = out_str + '\t' + str(rounded_val)
    out_f.write(out_str + '\n')

    # print(out_str+'\n')
    out_f.close()


def compute_path_based_effective_diffusion_from_R_M(R, W, alpha, top_pred_idx, pos_nodes_idx, alg_name):
    M = alg_alias[alg_name].get_M(W, alpha)
    R_diff_M = (R-M)
    effective_diff_all_targets_dict = {}
    for target in top_pred_idx:
        effective_diff_all_targets_dict[target] = sum(R_diff_M[target, pos_nodes_idx])/sum(R[target, pos_nodes_idx])

    return effective_diff_all_targets_dict

def compute_node_based_effective_diffusion_M_inv(M_inv, W, top_k_pred_idx, pos_nodes_idx,pred_scores):

    ######### the following code segement considers the M_inv matrix to be row normalized which it isn't.
    # Will not use this code for final version.
    # node_based_ed_all_targets_dict_test = {}
    # for target in top_k_pred_idx:
    #     target_score = sum(M_inv[pos_nodes_idx,target])
    #     neighbor_idx = set(np.nonzero(W.A.T[target])[0])
    #     non_nbr_pos = set(pos_nodes_idx) - neighbor_idx
    #     non_nbr_contr = sum(M_inv[list(non_nbr_pos), target])
    #     node_based_ed = non_nbr_contr / target_score
    #     node_based_ed_all_targets_dict_test[target] = node_based_ed
    #     print(target, '  ', node_based_ed)


    # node_based_ed_all_targets_dict_test={}
    # for target in top_k_pred_idx:
    #     target_score = sum(M_inv[target, pos_nodes_idx])
    #     neighbor_idx = set(np.nonzero(W.A[target])[0])
    #     non_nbr_pos = set(pos_nodes_idx) - neighbor_idx
    #     non_nbr_contr = sum(M_inv[target, list(non_nbr_pos)])
    #     node_based_ed = non_nbr_contr/target_score
    #     node_based_ed_all_targets_dict_test[target] = node_based_ed
    #     # print(target)

    ########### FINAL CODE
    node_based_ed_all_targets_dict = {}
    masked_M_inv = np.ma.masked_where(W.A>0,M_inv) #keep the non-neighbors only, so making the neighbors
    non_nbr_contr = masked_M_inv[:,pos_nodes_idx].sum(axis=1).data
    for target in top_k_pred_idx:
        node_based_ed_all_targets_dict[target] = non_nbr_contr[target]/pred_scores[target]
    ########################
    # for target in node_based_ed_all_targets_dict_test:
    #     d = abs(node_based_ed_all_targets_dict_test[target]-node_based_ed_all_targets_dict[target])
    #     assert d<10e-5,print('not same: ', d)

    return node_based_ed_all_targets_dict



def save_shortest_path_graphs(M_pathmtx,pos_nodes_idx, shortest_path_input_graph_file,
                              all_same_weight_input_graph_file):

    # taking transpose as newtorkx considers rows to be source and cols to be target but
    # in our case we normlaized M_pathmtx such a way that cols are sources and rows are target.
    G = nx.from_numpy_matrix(M_pathmtx.transpose(), create_using=nx.DiGraph())
    G1 = copy.deepcopy(G)

    # adding the supersource here
    ss_index = G.number_of_nodes()  # TODO make sure that in existng G I had max nodeidx = G.number_of_nodes - 1
    # graph_nodes = list(G.nodes)
    # graph_nodes.sort(reverse = True)
    # max_node_idx = graph_nodes[0]
    add_supersource(G, ss_index, pos_nodes_idx)

    # Now save this graph G as a list of edges with weights to pass it to java based k-sp algorithm
    nx.write_weighted_edgelist(G, shortest_path_input_graph_file)
    print('networkx graph creation done and saved in: ',shortest_path_input_graph_file )

    # now to find out actual number of shortest paths of less than a length l, we need the graph
    # G where all edges are of weight 1.
    nx.set_edge_attributes(G1, values=1, name='weight')
    # TODO add the supersource here
    add_supersource(G1, ss_index, pos_nodes_idx)
    nx.write_weighted_edgelist(G1, all_same_weight_input_graph_file)
    del G, G1

    print('networkx same weight graph creation done and saved in: ', all_same_weight_input_graph_file)
    return  ss_index

def main(config_map, k, **kwargs):
    """
    *config_map*: everything in the config file
    *kwargs*: all of the options passed into the script
    """

    # extract the general variables from the config map
    input_settings, input_dir, output_dir, alg_settings, kwargs \
        = config_utils.setup_config_variables(config_map, **kwargs)

    sig_cutoff = kwargs.get('stat_sig_cutoff')
    sig_str = "-sig%s" % (str(sig_cutoff).replace('.', '_')) if sig_cutoff else ""
    # m = kwargs.get('m')
    max_pathlen = kwargs.get('max_len')

    # for each dataset, extract the path(s) to the prediction files,
    # read in the predictions, and test for the statistical significance of overlap

    for dataset in input_settings['datasets']:
        print("Loading data for %s" % (dataset['net_version']))
        # load the network and the positive examples for each term
        net_obj, ann_obj, _ = setup_dataset(
            dataset, input_dir, **kwargs)
        prots, node2idx = net_obj.nodes, net_obj.node2idx

        for term in ann_obj.terms:
            term_idx = ann_obj.term2idx[term]
            orig_pos_idx, _ = alg_utils.get_term_pos_neg(ann_obj.ann_matrix, term_idx)
            orig_pos = [prots[p] for p in orig_pos_idx]
            pos_nodes_idx = [node2idx[n] for n in orig_pos if n in node2idx]
            n_pos = len(pos_nodes_idx)

            #If 'pos_k'=True, then the number of top predictions is equal to the number of positively annotated nodes
            # for this certain term.
            if kwargs.get('pos_k'):
                k = n_pos
                print('k: ', k)
            for alg_name in alg_settings:
                if (alg_settings[alg_name]['should_run'][0] == True) or (alg_name in kwargs.get('run_algs')):
                    # load the top predictions
                    print(alg_name)

                    if kwargs.get('balancing_alpha_only'): #in alg_setting[alg_name]['alpha'] put the balancing alpha
                        # get the balancing alpha for this network - alg - term
                        alpha_summary_filename = config_map['output_settings']['output_dir'] + \
                            "/viz/%s/%s/param_select/" % (dataset['net_version'], dataset[
                            'exp_name']) + '/' + alg_name + '/alpha_summary.tsv'
                        alpha_summary_df = pd.read_csv(alpha_summary_filename, sep='\t', index_col=None)[['term','balancing_alpha']]
                        term_2_balancing_alpha_dict = dict(zip(alpha_summary_df['term'], alpha_summary_df['balancing_alpha']))

                        balancing_alpha = term_2_balancing_alpha_dict[term]
                        alg_settings[alg_name]['alpha'] = [balancing_alpha]


                    alg_pred_files = config_utils.get_dataset_alg_prediction_files(
                        output_dir, dataset, alg_settings, [alg_name], **kwargs)
                    # get the alpha values to use
                    alphas = alg_settings[alg_name]['alpha']

                    for alpha, alg in zip(alphas, alg_pred_files):
                        pred_file = alg_pred_files[alg]
                        pred_file = script_utils.term_based_pred_file(pred_file, term)
                        if not os.path.isfile(pred_file):
                            print("Warning: %s not found. skipping" % (pred_file))
                            continue
                        print("reading %s for alpha=%s" % (pred_file, alpha))
                        df = pd.read_csv(pred_file, sep='\t')


                        # remove the original positives for downstream analysis
                        df = df[~df['prot'].isin(orig_pos)]
                        df.reset_index(inplace=True, drop=True)

                        if sig_cutoff:
                            df = config_utils.get_pvals_apply_cutoff(df, pred_file, **kwargs)

                        if k > len(df['prot']):
                            print("ERROR: k %s > num predictions %s. Quitting" % (k, len(df['prot'])))
                            sys.exit()

                        pred_scores = np.zeros(len(net_obj.nodes))
                        df = df[:k]
                        top_k_pred = df['prot']
                        top_k_pred_idx = [net_obj.node2idx[n] for n in top_k_pred]
                        pred_scores[top_k_pred_idx] = df['score'].values

                        # No need for including dataset['exp_name'] as the following matrix are seed node independent.
                        diff_mat_file = "%s/diffusion-mat-%s-a%s.npy" % (net_obj.out_pref, alg_name,
                                                                         str(alpha).replace('.', '_'))
                        fluid_flow_mat_file_M = "%s/fluid-flow-mat-%s-a%s.npy" % (net_obj.out_pref, alg_name,
                                                                                  str(alpha).replace('.', '_'))
                        fluid_flow_mat_file_R = "%s/fluid-flow-mat-R-%s-a%s.npy" % (net_obj.out_pref, alg_name,
                                                                                    str(alpha).replace('.', '_'))


                        ##########CREAE or LOAD THE DIFFUSION MAYTRICES
                        force_matrix = kwargs.get('force_matrix')
                        M_inv = alg_alias[alg_name].get_diffusion_matrix(net_obj.W, alpha=alpha,
                                diff_mat_file=diff_mat_file,force_run=force_matrix)
                        M_pathmtx, R = alg_alias[alg_name].get_fluid_flow_matrix(net_obj.W, alpha=alpha, \
                                        fluid_flow_mat_file_M=fluid_flow_mat_file_M, \
                                        fluid_flow_mat_file_R=fluid_flow_mat_file_R, force_run=force_matrix)

                        if alg_name =='rwr':
                            M_inv = M_inv * (alpha / float(len(pos_nodes_idx)))

                        ###COMPUTE NODE BASED EFFECTIVE DIFFUSION
                        target_to_node_based_ed = compute_node_based_effective_diffusion_M_inv(M_inv, net_obj.W,
                                                top_k_pred_idx, pos_nodes_idx,pred_scores)
                        print('node based ed done')

                        ###COMPUTE PATH BASED EFFECTIVE DIFFUSION USING MATRIX SUBTRACTION
                        target_to_path_based_ed = compute_path_based_effective_diffusion_from_R_M(R, net_obj.W,
                                                alpha, top_k_pred_idx, pos_nodes_idx, alg_name)
                        print('subtraction based path based ed done')

                        effective_diffusion_file = config_map['output_settings']['output_dir'] + "/viz/%s/%s/diffusion-analysis/%s/effective-diff-k%s-a%s%s.tsv" \
                                % (dataset['net_version'], term, alg_name, k, alpha, sig_str)
                        os.makedirs(os.path.dirname(effective_diffusion_file), exist_ok=True)


                        ed_df = pd.DataFrame()
                        assert list(target_to_path_based_ed.keys())==list(target_to_node_based_ed.keys()), print('list of targets doesnt match')

                        ed_df['prot_idx'] = list(target_to_path_based_ed.keys())
                        ed_df['prot'] = ed_df['prot_idx'].apply(lambda x: prots[x])
                        ed_df['path_based_ed'] = pd.Series(list(target_to_path_based_ed.values())).round(6)

                        ed_df['node_based_ed'] = pd.Series(list(target_to_node_based_ed.values())).round(6)

                        ed_df.to_csv(effective_diffusion_file, sep='\t', index=False)

                        print('Finished dataset: ', dataset, 'alg: ', alg_name)
                        del M_inv, R


if __name__ == "__main__":
    config_map, kwargs = parse_args()
    for k in kwargs.get('k_to_test'):
        main(config_map, k=k, **kwargs)
