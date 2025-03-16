#!/bin/sh -login
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=1G
#SBATCH --time=06:59:00
#SBATCH --job-name blastp.NCBI
#SBATCH --array 0-10



module load GCC/10.2.0  OpenMPI/4.0.5
module load BLAST+/2.11.0

# move the session to the working directory
cd ${SLURM_SUBMIT_DIR}/07_Annotation

# define split fasta file to use
if [[ ${SLURM_ARRAY_TASK_ID} -lt 10 ]]
then
    spl_file=protein_0${SLURM_ARRAY_TASK_ID}
else
    spl_file=protein_${SLURM_ARRAY_TASK_ID}
fi

# log split file names on
echo -en ${SLURM_ARRAY_TASK_ID} '\t' ${spl_file} '\n' >> spl_file_blast_log.txt


# blast search
blastp -db Drerio_GRCz11/ncbi_dataset/data/GCF_000002035.6/protein.faa -query Bbrach_ncbi_prots/${spl_file}.fasta -evalue 1e-10 -outfmt 6 -num_alignments 1 -soft_masking true -lcase_masking -max_hsps 1 -out split_blast_results/Bbrach.${spl_file}.NCBI.blastp



# monitor run and resources
# write job information to SLURM output file
scontrol show job $SLURM_JOB_ID 
# write resource usage to SLURM output file (powetools command)
js -j $SLURM_JOB_ID

