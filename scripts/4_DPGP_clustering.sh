# set variables
wd="/media/cdn-bc/RAID/Projects/FHyyy_Anya_ES/gitrepo/data"
env="py2.7"

##
cd $wd
conda activate $env

DP_GP_cluster.py -i matrices/vst_sig_cds_dmso_averaged.tsv \
-o DP_GP/vst_sig_cds --fast -n 1000 --max_iters 1000 \
--check_burnin_convergence --check_convergence \
--cluster_uncertainty_estimate --plot -p pdf

conda deactivate