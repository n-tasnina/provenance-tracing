#!/bin/bash
conda init
conda activate provenance

echo "compute scores for proteins by running algorithms (RL or RWR) with list of params mentioned in config file"
python src/FastSinkSource/run_eval_algs.py --num-pred-to-write -1 --config $1
echo "Select best parameter for RL and RWR that best balances two objectives: promoting label similarity among neighboring nodes and preserving alignment with each node’s initial label."
python src/scripts/prediction/alg_parameter_selection.py --config $1 
echo "compute node and path based effective diffusion"
python src/scripts/diffusion/effective_diffusion_node_path.py --force-contr --balancing-alpha-only  --pos-k --config $1
echo "Extract top contributing paths"
python src/scripts/diffusion/path_effective_diffusion_pathtype_supersink.py --balancing-alpha-only --pos-k --config $1
echo "Plot contribution coming along paths pf diffent lengths"
python src/scripts/diffusion/plot_path_type_based_diffusion.py --pos-k --balancing-alpha-only --config $1
echo "Plot node and path based effective diffusion across network and algorithms"
python src/scripts/diffusion/compare_across_networks_ed_types_balancing_alpha.py  --pos-k --config $1
echo "Compute diffusion betweenness score. Also compute the overlap between nodes with high betweenness score and viral interactors (and essential proteins)"
python src/scripts/betweenness/betweenness_src_spec.py --config $1
