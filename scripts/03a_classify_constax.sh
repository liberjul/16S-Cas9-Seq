#!/bin/bash --login
########## Define Resources Needed with SBATCH Lines ##########

#SBATCH --time=6:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --ntasks=1                 # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=20          # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem=32G                   # memory required per node - amount of memory (in bytes)
#SBATCH --job-name=16S-cas9-seq-constax        # you can give your job a name for easier identification (same as -J)
#SBATCH --output=%x-%j.SLURMout

########## Command Lines to Run ##########

cd ~/He_Lab/16S-Cas9-Seq/

for i in nd58 cas9d58 gp711
do
  constax -c 0.8 \
  -b -i ./exported-feature-table/rep-seqs_1_$i/dna-sequences.fasta \
  -n $SLURM_CPUS_PER_TASK \
  -d /mnt/ufs18/rs-022/bonito_lab/CONSTAX_May2020/SILVA_Bacteria_tf/SILVA_138_SSURef_tax_silva_bact.fasta \
  -f /mnt/ufs18/rs-022/bonito_lab/CONSTAX_May2020/SILVA_Bacteria_tf --mem $SLURM_MEM_PER_NODE -m 20 \
  -x ./taxonomy/tax_bact \
  -o ./taxonomy/out_bact \
  --high_level_db /mnt/ufs18/home-026/liberjul/Bonito_Lab/constax/test_new/SILVA138_NR_blastdb/SILVA138_NRdb.fasta \
  --conservative \
  --consistent \
  --msu_hpcc
  mv ./taxonomy/out_bact/constax_taxonomy.txt ./out_bact/constax_taxonomy_$i.txt
done

scontrol show job $SLURM_JOB_ID     ### write job information to output file
