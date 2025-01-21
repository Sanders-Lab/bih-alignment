#!/bin/bash
#
#SBATCH --job-name=aln_qc
#SBATCH --output=/data/cephfs-2/unmirrored/groups/sanders/projects/benedict/logs/aln_qc.txt
#
#SBATCH --cpus-per-task=64
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=1-00:00
#SBATCH --mem-per-cpu=2G

# script for running alignment QC *separately* to alignment pipeline (alignment must have already been run successfully)
# April 2024

echo "launching script at $(date)"
echo "SLURM_JOB_ID = $SLURM_JOB_ID, SLURM_JOB_NAME = $SLURM_JOB_NAME"

# INPUT: $1 = project name, $2  = organism
# Needs changing for each new project:
project_name=$1 # the name of the dir in /fast/groups/ag_sanders/work/data containig the reads (which should contain a dir named fastq/)
organism=$2

# Initiation
printf '\n ### Activating conda environment ####\n'

# dependencies: r samtools picard bwa fastqc multiqc bedtools
# the conda env can be installed with: conda env create -f //fast/groups/ag_sanders/work/tools/conda_envs/alignmentenv_20220505.yml

# test if conda env is present
condaenvexists=$(conda env list | grep 'alignmentenv' | awk '{print $1}')
[ -z ${condaenvexists} ] && { echo "ERROR: alignmentenv conda env not found in conda env list" ; exit ; } || echo ''
[ ${condaenvexists} = "alignmentenv" ] && echo "alignmentenv conda env found in conda env list" || { echo "ERROR: alignmentenv conda env not found in conda env list" ; exit ; }

conda init bash
source ~/.bashrc
conda activate alignmentenv

# test if sequencing project directory exists
if [[ ! -d ${SLURM_SUBMIT_DIR}/bam ]]
then
        echo "ERROR: this dir does not exist: ${SLURM_SUBMIT_DIR}/bam"
        echo "please launch from directory that contains a bam directory!"
        exit
fi

# test if slurm job was submitted from repo cloned from github
submit_dir_ls=$(ls $SLURM_SUBMIT_DIR | tr ' ' '\n' | grep 'bih-alignment' | wc -w)
if [[ ${submit_dir_ls} = 0 ]]
then
    echo "ERROR: the 'bih-alignment' directory could not be found in ${SLURM_SUBMIT_DIR}"
    echo "Please clone the github repo then submit the slurm job like this: 'sbatch bih-alignment/scripts/alignment_script.sh'"
    exit
fi
if [[ ! -d ${SLURM_SUBMIT_DIR}/bih-alignment/exec ]]
then
    echo "ERROR: the 'exec' directory could not be found in ${SLURM_SUBMIT_DIR}/bih-alignment/"
    echo "Do not alter the directory structure as it means the exec/ scripts cannot be found"
    exit
fi

# run script
echo "launching QC script ${SLURM_SUBMIT_DIR}/bih-alignment/exec/alignment_qc_exec.sh"
bash ${SLURM_SUBMIT_DIR}/bih-alignment/exec/alignment_qc_exec.sh \
	$project_name \
	$organism
