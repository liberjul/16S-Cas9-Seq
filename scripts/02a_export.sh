#!/bin/bash --login
########## Define Resources Needed with SBATCH Lines ##########

#SBATCH --time=6:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --ntasks=1                 # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=1          # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem=1G                   # memory required per node - amount of memory (in bytes)
#SBATCH --job-name=16S-cas9-seq_export_data        # you can give your job a name for easier identification (same as -J)
#SBATCH --output=%x-%j.SLURMout

########## Command Lines to Run ##########


cd ~/He_Lab/16S-Cas9-Seq/

conda activate qiime2-2020.8

for i in nd58 cas9d58 gp711
do
  # i=nd58
  # qiime metadata tabulate \
  #   --m-input-file feature_table/dada2-stats_1_$i.qza \
  #   --o-visualization feature_table/dada2-stats_1_$i.qzv
  #
  # # Visualize feature table 1 data summaries
  # qiime feature-table summarize \
  #   --i-table feature_table/table_1_$i.qza \
  #   --m-sample-metadata-file sample-metadata_$i.tsv \
  #   --o-visualization feature_table/table_1_$i.qzv
  #
  # qiime feature-table tabulate-seqs \
  #   --i-data feature_table/rep-seqs_1_$i.qza \
  #   --o-visualization feature_table/rep-seqs_1_$i.qzv

  qiime tools export \
    --input-path feature_table/rep-seqs_1_$i.qza \
    --output-path exported-feature-table/rep-seqs_1_$i

  qiime tools export \
    --input-path feature_table/table_1_$i.qza \
    --output-path exported-feature-table/table_1_$i
done

conda deactivate

scontrol show job $SLURM_JOB_ID     ### write job information to output file
