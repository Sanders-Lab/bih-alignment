#!/bin/bash
#
#SBATCH --job-name=aln_qc
#SBATCH --output=//fast/groups/ag_sanders/work/projects/benedict/logs/202206 .txt
#
#SBATCH --cpus-per-task=64
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --partition=highmem
#SBATCH --exclusive
#SBATCH --time=2-00:00
#SBATCH --mem-per-cpu=3G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=benedict.monteiro@mdc-berlin.de

# script for running alignment QC *separately* to alignment pipeline (alignment must have already been run successfully)
# June 2022

# Needs changing for each new project:
project_name= # the name of the dir in /fast/groups/ag_sanders/work/data containig the reads (which should contain a dir named fastq/)

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
if [[ ! -d /fast/groups/ag_sanders/work/data/${project_name} ]]
then
        echo "ERROR: this dir does not exist: /fast/groups/ag_sanders/work/data/${project_name}"
        echo "please set project_name in global options to a real directory!"
        exit
fi

# run script
echo "launching QC script //fast/groups/ag_sanders/work/projects/benedict/master_scripts/alignment/exec/alignment_qc_exec.sh"
bash //fast/groups/ag_sanders/work/projects/benedict/master_scripts/alignment/exec/alignment_qc_exec.sh \
	$project_name
