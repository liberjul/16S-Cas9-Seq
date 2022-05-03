#!/bin/bash --login
########## Define Resources Needed with SBATCH Lines ##########

#SBATCH --time=6:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --ntasks=1                 # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=1          # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem=1G                   # memory required per node - amount of memory (in bytes)
#SBATCH --job-name=16S-cas9-seq_filter_trim        # you can give your job a name for easier identification (same as -J)
#SBATCH --output=%x-%j.SLURMout

########## Command Lines to Run ##########

for i in nd58 cas9d58 gp711
do
  constax -c 0.8 -t \
  -b -i otus_R1_$i.fasta \
  -n $SLURM_CPUS_PER_TASK \
  -d /mnt/ufs18/rs-022/bonito_lab/CONSTAX_May2020/SILVA_Bacteria_tf/SILVA_138_SSURef_tax_silva_bact.fasta \
  -f /mnt/ufs18/rs-022/bonito_lab/CONSTAX_May2020/SILVA_Bacteria_tf --mem $SLURM_MEM_PER_NODE -m 20 \
  -x ./tax_bact \
  -o ./out_bact \
  --high_level_db /mnt/ufs18/home-026/liberjul/Bonito_Lab/constax/test_new/SILVA138_NR_blastdb/SILVA138_NRdb.fasta \
  --conservative \
  --consistent \
  --msu_hpcc
  mv ./out_bact/constax_taxonomy.txt ./out_bact/constax_taxonomy_$i.txt
done

scontrol show job $SLURM_JOB_ID     ### write job information to output file