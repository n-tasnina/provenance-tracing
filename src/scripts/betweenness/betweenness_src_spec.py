import copy
import sys
import numpy as np
import yaml
import argparse
import pandas as pd
import os

# sys.path.insert(1, "/data/tasnina/Provenance-Tracing/SARS-CoV-2-network-analysis/")
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../.."))
sys.path.insert(1, project_root)


# #use the following for running in terminal
from src.FastSinkSource.src.main import setup_dataset
from src.FastSinkSource.src.utils import config_utils
from src.FastSinkSource.src.algorithms import alg_utils
import src.scripts.betweenness.betweenness_utils as btns_utils
import src.scripts.betweenness.plot_utils  as btns_plot_utils
import src.scripts.utils  as script_utils
from scipy.sparse import eye, diags
from src.FastSinkSource.src.algorithms import rl_genemania as gm
# #
# from FastSinkSource.src.main import setup_dataset
# from FastSinkSource.src.utils import config_utils
# from FastSinkSource.src.algorithms import alg_utils
# import scripts.betweenness.betweenness_utils as btns_utils
# import scripts.betweenness.plot_utils  as btns_plot_utils
# import scripts.utils  as script_utils
from scipy.sparse import eye, diags
from src.FastSinkSource.src.algorithms import rl_genemania as gm

alg_plot_name = {'rwr': 'RWR', 'genemaniaplus': 'RL'}

def parse_args():
    parser = setup_opts()
    args = parser.parse_args()
    kwargs = vars(args)
    with open(args.config, 'r') as conf:
        config_map = yaml.load(conf, Loader=yaml.FullLoader)
    return config_map, kwargs


def setup_opts():
    ## Parse command line args.
    parser = argparse.ArgumentParser(description="Compute betweenness score for each node in the network")
    # general parameters
    group = parser.add_argument_group('Main Options')
    group.add_argument('--config', type=str, default="/home/grads/tasnina/Projects/Provenance-Tracing/fss_inputs/config_files/provenance/signor_s12.yaml"
                       , help="Configuration file used when running FSS. ")
    group.add_argument('--k-to-test', '-k', type=int, action='append', default=[332],
                       help="k-value(s) for which to get the top-k predictions to test. " +
                            "If not specified, will check the config file.")
    group.add_argument('--pos-k', action='store_true', default=True,
                       help="if true get the top-k predictions to test is equal to the number of positive annotations")

    group.add_argument('--essential-prot-file', type=str,
                       default ="fss_inputs/essential-prots/deg_annotation_e.csv",
                       help="This file should contain the essential genes for corresponding species")

    group.add_argument('--viral-prot-file', type=str,
                       default="fss_inputs/viral-prots/HVIDB.csv",
                       help="This file should contain the essential genes for corresponding species")
    # group.add_argument('--pleiotropic-prot-file', type=str,
    #                    default="/data/tasnina/Provenance-Tracing/SARS-CoV-2-network-analysis/"
    #                     "fss_inputs/pleiotropic-prots/pleiotropic_gene_human_ra.xls",
    #                    help="This file should contain the pleiotropic genes for corresponding species")
    # group.add_argument('--pleiotropic-corr', action='store_true', default=False,
    #                    help='If true we will find correlation between pleoitrophy and betweenness.')
    group.add_argument('--ks', type=str, action='append', default=[200,400,600,800,1000, 2000, 3000,4000, 5000,
                                                                   6000, 7000, 8000, 9000, 10000])
    group.add_argument('--run-algs', type=str, action='append', default=[])
    group.add_argument('--force-download', action='store_true', default=False,
                       help="Force re-downloading and parsing of the input files")
    group.add_argument('--force-run', action='store_true', default=False,
                       help="If true Compute betweenness even if it was computed before.")

    group.add_argument('--net-usage', action='store_true', default=True)
    #Nure: removing balancing_alpha_only as I will always pick balancing alpha
    # group.add_argument('--balancing-alpha-only', action='store_true', default=True,
    #                    help="Ignore alpha from config file rather take the alpha value\
    #                           that balanced the two loss terms in quad loss function for the corresponding\
    #                           network-term-alg")
    return parser


def main(config_map, **kwargs):


    """
    *config_map*: everything in the config file
    *kwargs*: all of the options passed into the script
    """
    # extract the general variables from the config map
    btns_ranks = kwargs.get('ks')
    force_run=kwargs.get('force_run')
    input_settings, input_dir, output_dir, alg_settings, kwargs \
        = config_utils.setup_config_variables(config_map, **kwargs)

    summary_file = output_dir + '/non_zero_contr_nodes.tsv'
    if kwargs.get('net_usage'):
        with open(summary_file, 'w') as f:
            f.write('alg_name\tdataset_name\ttotal_nodes\tnon_zero_contr_nodes\tfrac_non_zero_contr_nodes\n')
            f.close()

    #keep track of how many prots with non-zero betweenness appear as we consider contribution via
    #paths of len 2,3,4
    n_prots_appearing_at_each_pathlens={'network':[],'term':[],'alg':[],'alpha':[], 'pathlen_2':[],  'pathlen_3':[],  'pathlen_4':[]}

    #all_criteria_overlap_es_pvals_topks_multinet is a dict of dict of dict. first_key=alg_name, second_key = dataset/network name,
    # third_key=rank ks, value=overlap pvalue at each rank. The following dict will save overlap for org-specific essential genes and
    # sars2-- viral interactors set only.
    all_criteria_overlap_es_pvals_topks_multialg_multinet={alg_name:{} for alg_name in alg_settings}
    all_criteria_overlap_viral_pvals_topks_multialg_multinet={alg_name:{} for alg_name in alg_settings}

    #save the #ess and viral prots in the whole network after removing top_k preds and source proteins from list of prots
    ess_multinet={alg_name:{} for alg_name in alg_settings}
    viral_multinet={alg_name:{} for alg_name in alg_settings}

    dataset_names_str=''
    # non_src_non_top_prots_per_dataset = {}
    for dataset in input_settings['datasets']:
        print("Loading data for %s" % (dataset['net_version']))
        dataset_name = dataset['plot_exp_name']
        dataset_names_str = dataset_name if dataset_names_str =='' else dataset_names_str + '_' + dataset_name
        # load the network and the positive examples for each term
        net_obj, ann_obj, _ = setup_dataset(dataset, input_dir, **kwargs)
        prots, node2idx = net_obj.nodes, net_obj.node2idx

        total_net_nodes = len(prots)
        # get essential protein list.
        # map essential prots to uniprot ids
        ess_types = ['org', 'cell']
        viral_types = ['sars2--']

        ess_uniprots_dict = btns_utils.handle_essential_uniprot_mapping(**kwargs)
        #get human proteins that interact with viral prots
        viral_prot_file = os.path.join(project_root, kwargs.get('viral_prot_file'))
        viral_uniprots_dict =  btns_utils.parse_viral_prot_file(viral_prot_file)

        ##Directory for saving any betweenness related analysis result
        btns_out_dir = output_dir + '/betweenness/' + dataset['net_version'] + '/'

        for alg_name in alg_settings:
            if (alg_settings[alg_name]['should_run'][0] == True) or (alg_name in kwargs.get('run_algs')):
                print(alg_name)
                for term in ann_obj.terms:
                    alg_term_spec_btns_out_dir = btns_out_dir + alg_name + '/' + dataset['exp_name'] + '/'
                    term_idx = ann_obj.term2idx[term]
                    orig_pos_idx, _ = alg_utils.get_term_pos_neg(ann_obj.ann_matrix, term_idx)
                    orig_pos = [prots[p] for p in orig_pos_idx]
                    pos_nodes_idx = [node2idx[n] for n in orig_pos if n in node2idx]
                    assert len(orig_pos) == len(pos_nodes_idx), print('not all source present in net')
                    n_pos = len(pos_nodes_idx)

                    # If 'pos_k'=True, then the number of top predictions is equal to the number of positively annotated nodes
                    # for this certain term.
                    k=kwargs.get('k')
                    if kwargs.get('pos_k'):
                        k = n_pos
                        print('k: ', k)
                    # remove any positive protein present in viral interactors
                    # for viral_type in viral_uniprots_dict:
                    #     viral_uniprots_dict[viral_type] = viral_uniprots_dict[viral_type].difference(set(orig_pos))

                    # if kwargs.get('balancing_alpha_only'):
                    balancing_alpha = script_utils.get_balancing_alpha(config_map,dataset,alg_name,term)
                    alg_settings[alg_name]['alpha'] = [balancing_alpha]
                    #Get prediction file
                    alg_pred_files = config_utils.get_dataset_alg_prediction_files(
                        output_dir, dataset, alg_settings, [alg_name], **kwargs)
                    # get the alpha values to use
                    alphas = alg_settings[alg_name]['alpha']
                    count=0
                    for alpha, alg in zip(alphas, alg_pred_files):
                        count+=1

                        pred_file = alg_pred_files[alg]
                        pred_file = script_utils.term_based_pred_file(pred_file, term)
                        if not os.path.isfile(pred_file):
                            print("Warning: %s not found. skipping" % (pred_file))
                            continue
                        print("reading %s for alpha=%s" % (pred_file, alpha))
                        df = pd.read_csv(pred_file, sep='\t')
                        # # remove the original positives for downstream analysis
                        # df = df[~df['prot'].isin(orig_pos)]
                        # df.reset_index(inplace=True, drop=True)

                        if k > len(df['prot']):
                            print("ERROR: k %s > num predictions %s. Quitting" % (k, len(df['prot'])))
                            sys.exit()
                        df['prot_idx'] = df['prot'].apply(lambda x: net_obj.node2idx[x])


                        #COMPUTE BETWEENNESS
                        beta = alpha

                        #file to save the betweenness score and percent_rank(according to betweenness score ) of each gene
                        src_spec_btns_file = alg_term_spec_btns_out_dir +'btns_a' + str(alpha) + '.tsv'
                        #Compute BETWEENNESS score for each protein  and get a dataframe with proteins sorted in ascending order of btns score
                        sorted_src_spec_btns_df = btns_utils.handle_src_spec_btns\
                            (net_obj, alpha, prots,src_spec_btns_file, pos_nodes_idx, alg_name, force_run=force_run)

                        if kwargs.get('net_usage'):
                            #Analysis on impact of removing each non-src node from network
                            btns_matrix_file = src_spec_btns_file.replace('btns', 'btns_matrix').replace('.tsv','.pickle')
                            non_zero_contr_nodes = script_utils.global_net_usage(btns_matrix_file, df, k, orig_pos)

                            frac_non_zero_contr_nodes = non_zero_contr_nodes / total_net_nodes
                            with open(summary_file, 'a') as f:
                                f.write(f'{alg_name}\t{dataset_name}\t{total_net_nodes}\t{non_zero_contr_nodes}\t{frac_non_zero_contr_nodes}\n' )
                            f.close()



                        #compute and plot how many new prots appear with nonzero betweenness
                        # as we consider path lens of 2, 3, and 4. Did not filter out sources and top tragets yet

                        new_prots_appearing_at_each_pathlens = \
                            btns_utils.find_new_prots_appearing_at_each_pathlens(sorted_src_spec_btns_df)
                        n_prots_appearing_at_each_pathlens['network'].append(dataset_name)
                        n_prots_appearing_at_each_pathlens['term'].append(term)
                        n_prots_appearing_at_each_pathlens['alg'].append(alg_name)
                        n_prots_appearing_at_each_pathlens['alpha'].append(alpha)
                        n_prots_appearing_at_each_pathlens['pathlen_2'].append(copy.deepcopy(new_prots_appearing_at_each_pathlens['pathlen_2']))
                        n_prots_appearing_at_each_pathlens['pathlen_3'].append(copy.deepcopy(new_prots_appearing_at_each_pathlens['pathlen_3']))
                        n_prots_appearing_at_each_pathlens['pathlen_4'].append(copy.deepcopy(new_prots_appearing_at_each_pathlens['pathlen_4']))



                        #Now get the top predicted proteins
                        #here k= len(orig_pos)
                        pred_file = alg_pred_files[alg]
                        pred_file = script_utils.term_based_pred_file(pred_file, term)
                        top_k_predictions_df = script_utils.get_top_k_predictions(pred_file, alpha, len(orig_pos), orig_pos)

                        # filter the sorted_src_spec_btns_df to contain only the non-source, non-top predictions
                        sorted_filtered_src_spec_btns_df, sorted_df_pos, sorted_df_top_k = \
                            btns_utils.filter_sorted_src_spec_btns(sorted_src_spec_btns_df, top_k_predictions_df, orig_pos)


                        #include total number of non src, non top scoring proteins in the network to later get the random overlap between network nodes and essential/viral nodes.
                        btns_ranks_with_total_n = btns_ranks + [len(sorted_filtered_src_spec_btns_df['prot'].unique())]
                        ###Do analysis to find relationship between betweenness score and  ESSENTAIL PROT
                        # Kolmogorov-Smirnov test to see if essential proteins show significantly high
                        # btns score

                        # ****TODO uncomment later
                        ks_file = alg_term_spec_btns_out_dir + 'KS_pvals' + '_a' + str(alpha) + '.tsv'
                        KS_dict = btns_utils.handle_Kolmogorov_Smirnov_test(sorted_filtered_src_spec_btns_df, ess_uniprots_dict, viral_uniprots_dict)
                        script_utils.save_dict(KS_dict, ks_file)

                        marker= 'rank'
                        frac_prots_ge_btns_marker = btns_utils.prepare_plotdata_for_Kolmogorov_Smirnov\
                                                    (sorted_filtered_src_spec_btns_df, ess_uniprots_dict, viral_uniprots_dict,
                                                     marker =marker, ks=btns_ranks)

                        title = alg_plot_name[alg_name] + '_a_' + str(alpha) + '_' + term + '_' + dataset_name
                        ks_plt_file = alg_term_spec_btns_out_dir + 'KS_' +marker+ '_a' + str(alpha) + '.pdf'
                        # Plot for KS
                        btns_plot_utils.plot_KS(frac_prots_ge_btns_marker, marker, title, ks_plt_file)
                        # ****


                        # #*********************************  ESSENTIAL PROTS *****************************


                        print('\n\n ESSENTAIL PROTS ANALYSIS\n')
                        for ess_type in ess_types :
                                print(ess_type)
                                alg_term_ess_btns_corr_file = alg_term_spec_btns_out_dir + 'corr_ess_'+ess_type+'.tsv'
                                os.makedirs(os.path.dirname(alg_term_ess_btns_corr_file), exist_ok=True)

                                #Hypergeometric test
                                #compute overlap between top_ranked_prots and essential_prots. The ranking was done using
                                # paths of len 2,3,4 separately and all together.

                                all_criteria_overlap_es_pvals_topks = btns_utils.handle_Fishers_exact_test_in_topks\
                                                (sorted_filtered_src_spec_btns_df, ess_uniprots_dict[ess_type], btns_ranks_with_total_n)

                                #now compute  1. frac of src_nodes are essential 2. frac of predicted nodes by algorithms are essential
                                # 3. frac of nodes in networks excluding src and predicted_nodes are essential.
                                ess_in_pos = btns_utils.compute_frac_interesting_prot(sorted_df_pos['prot'], ess_uniprots_dict[ess_type])
                                ess_in_top = btns_utils.compute_frac_interesting_prot(sorted_df_top_k['prot'], ess_uniprots_dict[ess_type])
                                ess_in_net = btns_utils.compute_frac_interesting_prot(sorted_filtered_src_spec_btns_df['prot'], ess_uniprots_dict[ess_type])
                                if ess_type == 'org':
                                    all_criteria_overlap_es_pvals_topks_multialg_multinet[alg_name][dataset_name] = all_criteria_overlap_es_pvals_topks
                                    ess_multinet[alg_name][dataset_name]=ess_in_net

                                title = alg_plot_name[alg_name] + '_a_' + str(alpha) + '_' + term + '_' + dataset_name
                                overlap_pval_plt_file = alg_term_spec_btns_out_dir + 'overlap_ess_'+ess_type+'_a' + str(alpha) + '.pdf'
                                #Plot for hypergeometric test/Fisher's exact test
                                btns_plot_utils.plot_hypergeom_pval(all_criteria_overlap_es_pvals_topks, ess_in_pos,
                                            ess_in_top, ess_in_net, title, overlap_pval_plt_file, ks=btns_ranks)
                                btns_plot_utils.plot_hypergeom_pval(all_criteria_overlap_es_pvals_topks, ess_in_pos, ess_in_top,
                                ess_in_net, title,overlap_pval_plt_file.replace('.pdf', '_zoomed.pdf'), ks=[200,400,600,800,1000, 2000],
                                                                    rank_criteria=['betweenness'])
                                #Compute correlation between rank_percentiles of bins and percentage of essential prots in bins
                                pc_ess, pval_ess, mw_ess, prcntl_ess, prcnt_ess =\
                                    btns_utils.handle_percentile_percent_corr(sorted_filtered_src_spec_btns_df, ess_uniprots_dict[ess_type])
                                # SAVE correlations.
                                btns_utils.save_btns_corr(alg_term_ess_btns_corr_file, count, beta, pc_ess, pval_ess,mw_ess)


                                # Scatter plot for rank percentile and percentage of essential protein in each bin
                                ext_prcntl_ess, ext_prcnt_ess = btns_utils.find_interesting_prot_in_src_top(
                                    sorted_df_pos, sorted_df_top_k, ess_uniprots_dict[ess_type], prcntl_ess, prcnt_ess)
                                prcntl_prcnt_btns_ess_plt_file = alg_term_spec_btns_out_dir + 'scatter_ess_'+ess_type+'_a' + str(alpha) + '.pdf'
                                title = alg_plot_name[alg_name] + '_a_' + str(alpha) + '_' + term + '_' + dataset_name
                                x_label = 'percentile rank'
                                btns_plot_utils.scatter_plot(ext_prcntl_ess, ext_prcnt_ess, x_label=x_label,
                                                             y_label='percentage of essential prot : '+ess_type,
                                                             title=title, filename=prcntl_prcnt_btns_ess_plt_file)

                        # ****

                        # continue
                        #*********************************  VIRAL PROTS *****************************
                        print('\n\nVIRAL INTERACTOR ANALYSIS\n')
                        for viral_type in viral_types:
                            print(viral_type)

                            # Hypergeometric test
                            all_criteria_overlap_viral_pvals_topks = btns_utils.handle_Fishers_exact_test_in_topks \
                                (sorted_filtered_src_spec_btns_df, viral_uniprots_dict[viral_type], btns_ranks_with_total_n)

                            viral_in_pos = btns_utils.compute_frac_interesting_prot(sorted_df_pos['prot'],
                                                                                    viral_uniprots_dict[viral_type])
                            viral_in_top = btns_utils.compute_frac_interesting_prot(sorted_df_top_k['prot'],
                                                                                    viral_uniprots_dict[viral_type])
                            viral_in_net = btns_utils.compute_frac_interesting_prot(
                                sorted_filtered_src_spec_btns_df['prot'], viral_uniprots_dict[viral_type])

                            if viral_type == 'sars2--':
                                all_criteria_overlap_viral_pvals_topks_multialg_multinet[alg_name][dataset_name] = all_criteria_overlap_viral_pvals_topks
                                viral_multinet[alg_name][dataset_name] = viral_in_net
                            title = alg_plot_name[alg_name] + '_a_' + str(alpha) + '_' + term + '_' + dataset_name
                            overlap_pval_plt_file = alg_term_spec_btns_out_dir + 'overlap_viral_' + viral_type + '_a' + str(alpha) + '.pdf'

                            # Plot for hypergeometric test/Fisher's exact test
                            btns_plot_utils.plot_hypergeom_pval(all_criteria_overlap_viral_pvals_topks, viral_in_pos,
                                            viral_in_top, viral_in_net, title, overlap_pval_plt_file,ks=btns_ranks)

                            btns_plot_utils.plot_hypergeom_pval(all_criteria_overlap_viral_pvals_topks, viral_in_pos,
                                viral_in_top, viral_in_net, title, overlap_pval_plt_file.replace('.pdf', '_zoomed.pdf'),ks=[200,400,600,800,1000, 2000], rank_criteria=['betweenness'])

                            #Binwise Pearsons correlations
                            pc_viral, pval_viral, mw_viral, prcntl_viral, prcnt_viral = \
                                btns_utils.handle_percentile_percent_corr(sorted_filtered_src_spec_btns_df, viral_uniprots_dict[viral_type])
                            alg_term_viral_btns_corr_file = alg_term_spec_btns_out_dir + 'corr_viral_' + viral_type + '.tsv'

                            os.makedirs(os.path.dirname(alg_term_viral_btns_corr_file), exist_ok=True)
                            btns_utils.save_btns_corr(alg_term_viral_btns_corr_file, count, beta, pc_viral, pval_viral, mw_viral)

                            #Plot bin percentile-percentage scatter plot
                            ext_prcntl_viral, ext_prcnt_viral = btns_utils.find_interesting_prot_in_src_top(sorted_df_pos,
                                                                sorted_df_top_k, viral_uniprots_dict[viral_type], prcntl_viral, prcnt_viral)
                            # Plot a scatter plot of percentile-score and percentage of viral prot per bin
                            x_label = 'percentile rank'
                            prcntl_prcnt_btns_viral_plt_file = alg_term_spec_btns_out_dir + 'scatter_viral_'+viral_type+'_a' + str(alpha) + '.pdf'
                            btns_plot_utils.scatter_plot(ext_prcntl_viral, ext_prcnt_viral, x_label=x_label,
                                y_label='percentage of viral interactors : '+viral_type,
                                title=title, filename=prcntl_prcnt_btns_viral_plt_file)

                        #Plot binwise percentile rank for source proteins and top_k_predictions
                        src_top_rank_df = pd.DataFrame({'percent_rank': list(sorted_df_pos['percent_rank']) +
                            (list(sorted_df_top_k['percent_rank'])),'node_spec': (['pos'] * len(sorted_df_pos)) +
                            (['top_pred'] * len(sorted_df_top_k))})

                        src_top_rank_plot_file = alg_term_spec_btns_out_dir + 'percentile_rank_pos_top_a' + str(alpha) + '.pdf'
                        btns_plot_utils.box_plot(src_top_rank_df,x = 'node_spec', y='percent_rank', ymin=0, ymax=1,
                                 title=title, filename=src_top_rank_plot_file)

    btns_plot_utils.plot_prots_appearing_at_each_pathlens(n_prots_appearing_at_each_pathlens,
                                          filename= output_dir + '/betweenness/'+'new_appearing_prots.pdf')

    #Plot hypergeom pval plot for overlap between top_betwenness_prots and essential
    #**** TODO uncomment later
    overlap_es_pval_plt_file = output_dir + '/betweenness/' +dataset_names_str+ '_overlap_org_ess_'
    btns_plot_utils.lineplot_hypergeom_pval_multialg_multinet(all_criteria_overlap_es_pvals_topks_multialg_multinet,
                                                              ess_multinet,filename=overlap_es_pval_plt_file, ks=btns_ranks)
    btns_plot_utils.heatmap_hypergeom_pval_multialg_multinet(all_criteria_overlap_es_pvals_topks_multialg_multinet,
                                                             filename=overlap_es_pval_plt_file, )
    # ****

    overlap_viral_pval_plt_file = output_dir + '/betweenness/'+dataset_names_str + '_overlap_sars2--viral_'
    btns_plot_utils.lineplot_hypergeom_pval_multialg_multinet(all_criteria_overlap_viral_pvals_topks_multialg_multinet,
                                                              viral_multinet, filename=overlap_viral_pval_plt_file, ks=btns_ranks)
    btns_plot_utils.heatmap_hypergeom_pval_multialg_multinet(all_criteria_overlap_viral_pvals_topks_multialg_multinet,
                                                             filename=overlap_viral_pval_plt_file)

    btns_plot_utils.plot_prots_appearing_at_each_pathlens(n_prots_appearing_at_each_pathlens,
                                                          filename=output_dir + '/betweenness/' + 'new_appearing_prots.pdf')


    #plot bar chart for non_zero_contributing_nodes to top scoring nodes
    non_zero_contr_df = pd.read_csv(summary_file, sep='\t')
    btns_plot_utils.plot_non_zero_contr_nodes(non_zero_contr_df, output_dir)




if __name__ == "__main__":
    config_map, kwargs = parse_args()
    main(config_map, **kwargs)
