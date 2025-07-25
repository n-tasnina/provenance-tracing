import os, sys
import yaml
import argparse
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from datetime import datetime

sys.path.insert(1, "/data/tasnina/Provenance-Tracing/SARS-CoV-2-network-analysis/")
from src.FastSinkSource.src.main import setup_dataset
from src.FastSinkSource.src.utils import config_utils
from src.FastSinkSource.src.algorithms import alg_utils
import src.scripts.utils as script_utils
from src.scripts.plot_utils import *
import json
import copy

from src.graphspace import post_to_graphspace_base as gs
from src.graphspace.signor_gs_utils import *
from src.graphspace import backend_utils as back_utils
from src.graphspace import post_to_graphspace_wrapper as wrapper
from src.utils import network_processing_utils as net_utils
from src.utils import mapping_utils as map_utils


pathway_alias = {
    'SIGNOR-SC-Apoptosis_09_07_23.tsv' : 'Apoptosis',
    'SIGNOR-SC-Attachment-and-Entry_09_07_23.tsv': 'Attachment-and-Entry',
    'SIGNOR-SC-Cytokyne-Storm_09_07_23.tsv': 'Cytokyne-Storm',
    'SIGNOR-SC-Inflammatory-Response_09_07_23.tsv': 'Inflammatory-Response',
    'SIGNOR-SC-MAPK-Pathway_09_07_23.tsv': 'MAPK-Pathway',
    'SIGNOR-SC-Stress-Granules_09_07_23.tsv': 'Stress-Granules',
    'SIGNOR-SCFI_09_07_23.tsv': 'Fibrosis',
    'SIGNOR-SCISOVG_09_07_23.tsv': 'Innate Response to dsRNA',
    'SIGNOR-SCUP_09_07_23.tsv': 'ER Stress'
}

def parse_args():
    parser = setup_opts()
    args = parser.parse_args()
    kwargs = vars(args)
    with open(args.config, 'r') as conf:
        config_map = yaml.load(conf, Loader=yaml.FullLoader)
        # config_map = yaml.load(conf)
    return config_map, kwargs


def setup_opts():
    ## Parse command line args.
    parser = argparse.ArgumentParser(description="Script to visualize top contributing paths to k top scoring predictions")
    # general parameters
    group = parser.add_argument_group('Main Options')

    #ALGO specific arguments
    group.add_argument('--config', type=str, default="/data/tasnina/Provenance-Tracing/SARS-CoV-2-network-analysis/"
                        "fss_inputs/config_files/provenance/signor_s12.yaml"
                       , help="Configuration file used when running FSS. ")
    group.add_argument('--id-mapping-file', type=str, default = "/data/tasnina/Provenance-Tracing/SARS-CoV-2-network-analysis/"
                        "datasets/mappings/human/uniprot-reviewed-status.tab.gz",
                       help="Table downloaded from UniProt to map to gene names. Expected columns: 'Entry', 'Gene names', 'Protein names'")
    group.add_argument('--evidence-file', default = "/data/tasnina/Provenance-Tracing/SARS-CoV-2-network-analysis/"
                        "datasets/networks/signor-cc/all_data_22_12_22.tsv", help='File containing the evidence for each edge in the interactome.')
    group.add_argument('--k-to-test', '-k', type=int, action='append', default=[332],
                       help="k-value(s) for which to get the top-k predictions to test. " +
                            "If not specified, will check the config file.")
    group.add_argument('--pos-k', action='store_true', default=True,
                       help="if true get the top-k predictions to test is equal to the number of positive annotations")
    group.add_argument('--n-sp', '-n', type=int, default=1000,
                       help="n-sp is the number of shortest paths to be considered" +
                            "If not specified, will check the config file. Default=20")
    group.add_argument('--n-sp-viz', type=int, default=100,
                       help="How many top paths to vizualize" +
                            "Default=20")

    group.add_argument('--balancing-alpha-only', action='store_true', default=True,
                       help="Ignore alpha from config file rather take the alpha value\
                            that balanced the two loss terms in quad loss function for the corresponding\
                            network-term-alg")


    ############### GRAPH ATTR

    group.add_argument('--force-attr', action='store_true', default=True,
                       help="If true create new graph_attr_file.")


    #GRAPHSPACE POSTING

    group.add_argument('-U', '--username', type=str,  default='tasnina@vt.edu',
                      help='GraphSpace account username to post graph to. Required')
    group.add_argument('-P', '--password', type=str, default='1993Hello#GraphSpace',
                      help='Username\'s GraphSpace account password. Required')
    group.add_argument('--graph-name', type=str,  default='signor',
                      help='Graph name for posting to GraphSpace. Default: "test".')
    # replace with an option to write the JSON file before/after posting
    # parser.add_option('', '--outprefix', type='string', metavar='STR', default='test',
    #                  help='Prefix of name to place output files. Required.')
    group.add_argument('--group', type=str,
                      help='Name of group to share the graph with.')
    group.add_argument('--make-public', action="store_true", default=False,
                      help='Option to make the uploaded graph public.')
    # TODO implement and test this option
    # parser.add_argument('', '--group-id', type='string', metavar='STR',
    #                  help='ID of the group. Could be useful to share a graph with a group that is not owned by the person posting')
    group.add_argument( '--tags', type=str,action="append",
                      help='Tag to put on the graph. Can list multiple tags (for example --tag tag1 --tag tag2)')
    group.add_argument( '--apply-layout', type=str,
                      help='Specify the name of a graph from which to apply a layout. Layout name specified by the --layout-name option. ' +
                           'If left blank and the graph is being updated, it will attempt to apply the --layout-name layout.')
    group.add_argument( '--layout-name', type=str, default='layout1',
                      help="Name of the layout (of the graph specified by the --apply-layout option). " +
                           "X and y coordinates of nodes from that layout will be applied to matching node IDs in this graph. Default: 'layout1'")
    group.add_argument('--parent-nodes', action="store_true", default=True,
                      help='Use parent/group/compound nodes for the different node types')

    group.add_argument('--out-pref', type=str, metavar='STR',default='/data/tasnina/Provenance-Tracing/SARS-CoV-2-network-analysis/outputs/'
                                                                     'graphspace/',
                      help='Prefix of name to place output files. ')
    return parser



def get_edge_to_rank_map(paths_df):
    '''
    return a dict with key = tuple indicating (a,b) edge, valu = list containing rank of each path(in top n_sp
    (e.g., 1000) paths ) where this edges appeared
    '''
    edge_to_rank_mapping = {}
    for rank, row in paths_df.iterrows():
        edges_in_path = net_utils.path_nodes_to_set_of_edges(row['path_prots'])
        for edge in edges_in_path:
            if not (edge in edge_to_rank_mapping):
                edge_to_rank_mapping[edge] = [rank+1]
            else:
                edge_to_rank_mapping[edge].append(rank+1)
    return edge_to_rank_mapping

def process_causal_pathways(pathway_file):
    '''
    return: A set of tuples where each tuple (a,b) is an edge from a to b.
    '''
    cov_pathways = pd.read_csv(pathway_file, sep='\t')
    #todo: the following filtered out half the interaction. need to make a decision about
    #what is the right way to preprocess.
    cov_pathways = cov_pathways[(cov_pathways['TYPEA'] == 'protein')
                                & (cov_pathways['TYPEB'] == 'protein')
                                &(cov_pathways['TAX_ID'] == 9606)]
    cov_pathways = cov_pathways[['IDA', 'IDB']]
    unique_proteins = set(cov_pathways['IDA']).union(set(cov_pathways['IDB']))
    print('number of proteins in processed COVID causal pathways: ', len(unique_proteins))

    causal_edges = set(zip(cov_pathways['IDA'], cov_pathways['IDB']))
    print('number of edges in processed COVID causal pathways: ', len(causal_edges))

    return causal_edges



def process_top_paths_with_frac_causal(frac_causal_in_each_path_df):
    '''
    Function: Extract information on 'interesting' top contributing paths. And assign each path to a
    hallmark SARS-CoV causal pathway.
    Return: A dict where key=name of hallmark pathway, value = list of list where inner list is the nodes along a path
    that has been assigned to the corresponding hallmark pathway.
    '''
    #filter out paths of length 1
    frac_causal_in_each_path_df = frac_causal_in_each_path_df[frac_causal_in_each_path_df['length']>1]

    #keep interesting paths only, i.e., a certain fraction of edges is from COVID causal pathway. For now I am considering
    #paths with 2/2, 3/3, 2/3, 4/4/, 3/4, 5/5, 4/5 where numerator = #overlapping edges between COVID-causal and path, denominator=#edges in the path.
    frac_causal_in_each_path_df = frac_causal_in_each_path_df[frac_causal_in_each_path_df['COVID']>0.65]

    #Assign each remaining path to the hallmark SARS-CoV pathway with which it has the highest overlap. If multiple pathways
    # have the same overlap value then assign the path to all of them.

    hallmark_pathways = list(frac_causal_in_each_path_df.columns[6:])
    hallmark_pathways_path_assignment_dict = {key: [] for key in hallmark_pathways}
    for index, row in frac_causal_in_each_path_df.iterrows():
        max_overlap = 0
        assigned_pathway = []
        #find out the highest overlapping pathway(s) for a path
        for pathway in hallmark_pathways_path_assignment_dict:
            if row[pathway]>max_overlap:
                assigned_pathway = [pathway]
                max_overlap = row[pathway]
            elif row[pathway]==max_overlap:
                assigned_pathway.append(pathway)

        #save the path under the highest overlapping pathway(s)
        for pathway in assigned_pathway:
            hallmark_pathways_path_assignment_dict[pathway].append(row['path'])
    return hallmark_pathways_path_assignment_dict

def main(config_map, k, **kwargs):
    """
    *config_map*: everything in the config file
    *kwargs*: all of the options passed into the script
    """
    nsp = kwargs.get('n_sp')
    sig_cutoff = kwargs.get('stat_sig_cutoff')
    sig_str = "-sig%s" % (str(sig_cutoff).replace('.', '_')) if sig_cutoff else ""

    input_settings, input_dir, output_dir, alg_settings, kwargs \
        = config_utils.setup_config_variables(config_map, **kwargs)
    # gene name for corresponding uniprot is the node label
    if kwargs.get('id_mapping_file') is not None:
        uniprot_to_gene, uniprot_to_protein_names = map_utils.load_gene_names(kwargs.get('id_mapping_file'))

    for dataset in input_settings['datasets']:
        net_obj, ann_obj, _ = setup_dataset(dataset, input_dir, **kwargs)
        prots, node2idx = net_obj.nodes, net_obj.node2idx
        prot_universe = set(prots) #TODO discuss with Murali about what prot set to use as universe
        dataset_name = config_utils.get_dataset_name(dataset)

        print("\t%d prots in universe" % (len(prot_universe)))
        for term in ann_obj.terms:
            term_idx = ann_obj.term2idx[term]
            orig_pos_idx, _ = alg_utils.get_term_pos_neg(ann_obj.ann_matrix, term_idx)
            orig_pos = [prots[p] for p in orig_pos_idx]
            targets = list(set(prots).difference(set(orig_pos))) #target = nodes in network that are not EDNs

            pos_nodes_idx = [node2idx[n] for n in orig_pos if n in node2idx]
            n_pos = len(pos_nodes_idx)
            # If 'pos_k'=True, then the number of top predictions is equal to the number of positively annotated nodes
            # for this certain term.
            if kwargs.get('pos_k'):
                k = n_pos
                print('k: ', k)
            for alg_name in alg_settings:
                if (alg_settings[alg_name]['should_run'][0] == True):
                    # load the top predictions
                    print(alg_name)
                    if kwargs.get('balancing_alpha_only'):
                        balancing_alpha = script_utils.get_balancing_alpha(config_map, dataset, alg_name, term)
                        alg_settings[alg_name]['alpha'] = [balancing_alpha]
                    alphas = alg_settings[alg_name]['alpha']
                    alg_pred_files = config_utils.get_dataset_alg_prediction_files(
                        output_dir, dataset, alg_settings, [alg_name], **kwargs)

                    for alpha, alg in zip(alphas, alg_pred_files):

                        #Read file with top contributing paths
                        nsp_processed_paths_file = config_map['output_settings'][
                                                    'output_dir'] + "/viz/%s/%s/diffusion-path-analysis/%s/shortest_path_2ss/" \
                                                    "processed_shortest-paths-2ss-k%s-nsp%s-a%s%s.tsv" % (
                                                    dataset['net_version'], term, alg_name,k, nsp, alpha, sig_str)
                        # reading 'path_prots' column value as list
                        paths_df = pd.read_csv(nsp_processed_paths_file, sep='\t', index_col=None,
                                               converters={'path_prots': pd.eval})
                        paths = list(paths_df['path_prots'])  # this is one nested list. flatten it.
                        edge_to_rank_map = get_edge_to_rank_map(paths_df)


                        # pred score file
                        pred_file = alg_pred_files[alg]
                        pred_file = script_utils.term_based_pred_file(pred_file, term)
                        df = pd.read_csv(pred_file, sep='\t')
                        df_nonpos = df[~df['prot'].isin(orig_pos)].reset_index(drop=True)
                        # top_targets = nodes that have high predicted score by diffusion alg
                        top_targets = list(df_nonpos['prot'])[0:k]

                        kwargs['graph_name'] = dataset_name +'-top-paths-' + \
                            str(nsp) + '-causal_paths-'+ str(datetime.now())

                        # ************************** get edges in COVID-19 causal network
                        causal_pathways_dict = {}
                        pathway_dir = '/data/tasnina/Provenance-Tracing/SARS-CoV-2-network-analysis' \
                                      '/datasets/networks/signor-cc/COVID-causal/'
                        pathway_file = pathway_dir + 'SIGNOR-CCN_03_07_23.tsv'
                        causal_pathways_dict['COVID'] = process_causal_pathways(pathway_file)

                        # get info on SARS-CoV Hallmark pathways
                        file_dir_list = os.listdir(pathway_dir)
                        for f in file_dir_list:
                            if '-SC' in f:  # all SARS-CoV pathway file have '-SC' substring
                                causal_pathways_dict[pathway_alias[f]] = process_causal_pathways(pathway_dir + f)

                        # Read information on association between top contributing paths and SARS-COV causal paths
                        top_paths_causal_file = output_dir + '/causal_pathway/each_top_path_and_causal_pathways.tsv'
                        frac_causal_in_each_path_df = pd.read_csv(top_paths_causal_file, sep='\t',
                                                                  index_col=None, converters={'path': pd.eval})
                        hallmark_pathways_path_assignment_dict = process_top_paths_with_frac_causal(
                            copy.deepcopy(frac_causal_in_each_path_df))

                        for pathway in hallmark_pathways_path_assignment_dict:
                            if len(hallmark_pathways_path_assignment_dict[pathway])>0:
                                assigned_path_edges = net_utils.get_edges_from_list_of_paths(hallmark_pathways_path_assignment_dict[pathway])
                                assigned_path_edges_rank_dict = {}
                                for edge in assigned_path_edges:
                                    assigned_path_edges_rank_dict[edge] = edge_to_rank_map[edge]

                                pathway_edges = causal_pathways_dict[pathway]
                                kwargs['graph_name'] = 'causal_pathway:'+pathway+'-'+str(datetime.now())
                                wrapper.call_post_to_graphspace(assigned_path_edges_rank_dict, orig_pos, targets, top_targets,
                                    uniprot_to_gene, uniprot_to_protein_names,
                                    output_dir + '/graphspace/color.txt', causal_edges = pathway_edges, **kwargs)


if __name__ == "__main__":
    config_map, kwargs = parse_args()
    for k in kwargs.get('k_to_test'):
        main(config_map, k=k, **kwargs)

