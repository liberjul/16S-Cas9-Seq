#!/bin/bash --login
########## Define Resources Needed with SBATCH Lines ##########

#SBATCH --time=01:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --ntasks=1                 # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=1          # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem=10G                   # memory required per node - amount of memory (in bytes)
#SBATCH --job-name=16S-cas9-seq_filter_trim_cluster        # you can give your job a name for easier identification (same as -J)
#SBATCH --output=%x-%j.SLURMout

########## Command Lines to Run ##########

cd ~/He_Lab/16S-Cas9-Seq/

conda activate qiime2-2020.8

for i in nd58 cas9d58 gp711
do
  qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path ./data/raw_fastq/$i \
    --input-format CasavaOneEightSingleLanePerSampleDirFmt \
    --output-path demux-paired-end-$i.qza

  qiime dada2 denoise-paired \
    --i-demultiplexed-seqs demux-paired-end-$i.qza \
    --p-trunc-len-f 235\
    --p-trunc-len-r 235\
    --p-trim-left-f 19 \
    --p-trim-left-r 18 \
    --p-n-threads 40 \
    --verbose \
    --o-denoising-stats feature_table/dada2-stats_1_$i.qza \
    --o-representative-sequences feature_table/rep-seqs_1_$i.qza \
    --o-table feature_table/table_1_$i.qza
done

conda deactivate

scontrol show job $SLURM_JOB_ID     ### write job information to output file
