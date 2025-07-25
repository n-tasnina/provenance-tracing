#!/bin/bash
conda activate sarscov2-net-new

python src/FastSinkSource/run_eval_algs.py --num-pred-to-write -1 --config $1
python src/scripts/prediction/alg_parameter_selection.py --config $1
python src/scripts/diffusion/effective_diffusion_node_path.py --force-contr --balancing-alpha-only  --pos-k --config $1
python src/scripts/diffusion/path_effective_diffusion_pathtype_supersink.py --balancing-alpha-only --pos-k --config $1
python src/scripts/diffusion/plot_path_type_based_diffusion.py --pos-k --balancing-alpha-only --config $1
python src/scripts/diffusion/compare_across_networks_ed_types_balancing_alpha.py  --pos-k --config $1
python src/scripts/betweenness/betweenness_src_spec.py --config $1 --master-config $2
