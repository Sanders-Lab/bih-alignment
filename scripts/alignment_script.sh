#!/bin/bash
#
#SBATCH --cpus-per-task=80
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=2-00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --exclusive
#SBATCH --partition=highmem

# BIH cluster Paired end sequencing data alignment script July 2023
# This script takes .fastq format files, performs QC on them, before aligning to reference genome and outputting alignment QC metrics 
echo 'Running alignment script' ; date

##################################################################################################
# 1. Set global options
##################################################################################################
printf '\n ### 1. Setting global options ####\n'

# INPUT: $1 = project name [P3028], $2  = mate1_suffix [.1.fastq.gz], $3 = organism [human/mouse]
# should be submitted from a directory containing a directory of fastq files of interest

# Needs changing these for each new project:
project_name=$1 # the name of the Sample
mate1_suffix=$2 # suffix of read pair mate 1 fastq file (e.g. .1.fastq.gz)
organism=$3 # whether alignment is on human or mouse data
memperthread=2G # memory assigned per thread in slurm job
run_qc=TRUE # whether the alignment QC script should be run automatically [TRUE/FALSE]

# Probably don't need changing at the start of a new project:
mate2_suffix=$(echo $mate1_suffix | sed 's/1/2/') # suffix of mate 2 fastq file
fastq_dir=$SLURM_SUBMIT_DIR/fastq # dir containing fastq format files
n_threads=$(nproc) # number of threads given to slurm job, for highest efficiency use a multiple of 4 (e.g. 32, 64, etc.)
n_threads_divided=$(expr $n_threads / 4)
submit_dir=$SLURM_SUBMIT_DIR

# set reference genome
if [ $organism = 'human' ]
then
	ref_genome=/fast/groups/ag_sanders/work/data/references/genomes/human/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna # path to reference genome
	echo "Aligning human data, reference genome set to: $ref_genome"
elif [ $organism = 'human_t2t' ]
then
	ref_genome=/fast/groups/ag_sanders/work/data/references/genomes/human/t2t/hs1.fa.gz # path to reference genome
	echo "Aligning human T2T data, reference genome set to: $ref_genome"
 	organism=human # change for downstream to treat T2T same as hg38 human
elif [ $organism = 'mouse' ]
then
	ref_genome=/fast/groups/ag_sanders/work/data/references/genomes/mouse_mm39/mm39.fa.gz # path to reference genome
	echo "Aligning mouse data, reference genome set to: $ref_genome"
else
	echo "ERROR: command line argument 3 must be one of [human/human_t2t/mouse], currently it is: $organism"
	exit
fi


# print to log
echo "aligning fastq files in $fastq_dir"
echo "project_name: ${project_name}, mate1_suffix: ${mate1_suffix}, mate2_suffix: ${mate2_suffix}"

##################################################################################################
# 2. Activate conda environment
##################################################################################################
printf '\n ### 2. Activating conda environment ####\n'

# dependencies: r samtools picard bwa fastqc multiqc bedtools
# the conda env can be installed with: conda env create -f //fast/groups/ag_sanders/work/tools/conda_envs/alignmentenv_20220505.yml

# test if conda env is present
condaenvexists=$(conda env list | grep 'alignmentenv' | awk '{print $1}')
[ -z ${condaenvexists} ] && { echo "ERROR: alignmentenv conda env not found in conda env list" ; echo "please install the environment following the instructions in the GitHub README" ; exit ; } || echo 'conda environment found successfully, activating now'
# [ ${condaenvexists} = "alignmentenv" ] && echo "alignmentenv conda env found in conda env list" || { echo "ERROR: alignmentenv conda env not found in conda env list" ; exit ; }

conda init bash
source ~/.bashrc
conda activate alignmentenv

##################################################################################################
# 3. Initiation
##################################################################################################
printf '\n ### 3. Initiating script #####\n'

# test if sequencing project directory exists
if [[ ! -d $SLURM_SUBMIT_DIR ]]
then
	echo "ERROR: this dir does not exist: $SLURM_SUBMIT_DIR"
	echo "please set command line option 1 as a real directory!"
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

# confirm that there are .fastq files in the fastq_dir
[ ! $(ls ${fastq_dir}/*${mate1_suffix} | wc -l) -ge 1 ] && { echo "ERROR: no files were found with the suffix ${mate1_suffix} in ${fastq_dir}" ; echo "Please change the mate1_suffix manually in ${SLURM_SUBMIT_DIR}/scripts/alignment_script.sh" ; exit ; }
testfile=$(ls $fastq_dir | head -n1)
if [ ! $(zcat ${fastq_dir}/${testfile} | cut -c 1 | head -n1) = '@' ]
then
	echo "Error: No fastq files were found in ${fastq_dir}/, exiting script"
    echo "(the first character of the file ${fastq_dir}/${testfile} is not '@', suggesting it is not correct fastq format)"
	exit
fi

# create directories
tmp_dir=/fast/groups/ag_sanders/scratch/sequencing_tmp/${project_name} ; mkdir -p -m 775 $tmp_dir
bam_dir=${SLURM_SUBMIT_DIR}/bam; mkdir -m 775 $bam_dir
qc_dir=${SLURM_SUBMIT_DIR}/qc ; mkdir -m 775 $qc_dir
statsdir=${qc_dir}/alignment_stats ; mkdir -m 775 $statsdir
logsdir=${SLURM_SUBMIT_DIR}/logs ; mkdir -m 775 $logsdir

echo "Temporary/intermediate files will be written to ${tmp_dir}"
echo "Final .bam files will be written to ${bam_dir}"



# test if $run_qc is set correctly
if [ $run_qc = 'TRUE' ]
then
	echo "QC script will be run following alignment completion"
elif [ $run_qc = 'FALSE' ]
then
	echo "QC script will *NOT* be run following alignment completion"
else
	echo "ERROR: run_qc in Global Options must = TRUE or FALSE (currently it is neither!)"
	exit
fi

##################################################################################################
# 4. Run FastQC quality assessment on raw sequencing data
##################################################################################################
printf '\n ### 4. Run QC on sequencing reads #####\n'

fastqc_dir=${qc_dir}/fastqc ; mkdir -m 775 $fastqc_dir
multiqc_dir=${qc_dir}/multiqc ; mkdir -m 775 $multiqc_dir

# run FastQC to generate quality metrics on
echo "Running FastQC"
fastqc -q -t $n_threads -o ${fastqc_dir} ${fastq_dir}/*fastq.gz

# run multiqc
echo "Running MutiQC to combine FastQC results"
multiqc -o $multiqc_dir $fastqc_dir
echo "MutiQC results saved to: ${multiqc_dir}"


##################################################################################################
# 5. Alignment
##################################################################################################
printf '\n ### 5. Running alignment #####\n'

# make unique list of library names by removing suffixes
libraries=$(ls $fastq_dir/*fastq.gz \
	| sed -e "s?${fastq_dir}/??g;s?${mate1_suffix}??g;s?${mate2_suffix}??g" \
	| sort -u)

samdir=${tmp_dir}/sam ; mkdir -m 775 $samdir

alnlog=${logsdir}/$(date +%Y%m%d)_${project_name}_bwamem_log.txt
echo "see ${alnlog} for detailed BWA MEM log"

# perform alignment all libraries
for library in $libraries
do
	(
	echo "performing alignment on ${library}"
	echo "performing alignment on ${library}" >> $alnlog
	input1=${fastq_dir}/${library}${mate1_suffix}
	input2=${fastq_dir}/${library}${mate2_suffix}

	bwa mem -t 8 -v 3 -R $(echo "@RG\tID:${library}\tSM:${project_name}") \
    		$ref_genome $input1 $input2 > ${samdir}/${library}.sam 2>> $alnlog
	) &
	if [[ $(jobs -r -p | wc -l) -ge 10 ]]; # allows n_threads / 4 number of iterations to be executed in parallel
	then
		wait -n # if there are n_threads_divided iterations running wait here for space to start next job
	fi
done
wait # wait for all jobs in the above loop to be done

echo "Alignment completed, SAM files are saved to: ${samdir}"


##################################################################################################
# 6. Process BAM files
##################################################################################################
printf '\n ### 6. Processing BAM files #####\n'

# make dir on scratch drive for intermediate bam files
bam_tmpdir=${tmp_dir}/bam_tmpdir ; mkdir -m 775 $bam_tmpdir
mdup_metrics_dir=${statsdir}/mdup_metrics ; mkdir -m 775 $mdup_metrics_dir
mdup_tmpdir=${tmp_dir}/mdup_tmp ; mkdir -m 775 $mdup_tmpdir

duplog=${logsdir}/$(date +%Y%m%d)_${project_name}_mdup_log.txt
echo "see ${duplog} for detailed BWA MEM log"

# convert, process, filter SAM files
for library in $libraries
do
    (
	echo "processing ${library}.sam"
	# convert SAM to BAM
	samtools view -@ 3 -h -b ${samdir}/${library}.sam > ${bam_tmpdir}/${library}.bam # convert SAM to BAM

	# sort BAM file
	samtools sort -@ 3 -m $memperthread ${bam_tmpdir}/${library}.bam > ${bam_tmpdir}/${library}.sort.bam
	samtools index -@ 3 ${bam_tmpdir}/${library}.sort.bam # generate index

	# Mark duplicated reads in BAM file
	picard MarkDuplicates -I ${bam_tmpdir}/${library}.sort.bam \
        	-O ${bam_tmpdir}/${library}.sort.mdup.bam \
		-M ${mdup_metrics_dir}/${library}_mdup_metrics.txt \
        	--QUIET true --VERBOSITY ERROR --TMP_DIR ${mdup_tmpdir} 2>> $duplog
	samtools index -@ 3 ${bam_tmpdir}/${library}.sort.mdup.bam # generate index

	# copy sorted marked duplicates BAM files and their indexes to work drive
	cp ${bam_tmpdir}/${library}.sort.mdup.bam* $bam_dir
	) &
	if [[ $(jobs -r -p | wc -l) -ge $n_threads_divided ]]; # allows n_threads / 4 number of iterations to be executed in parallel
    	then
		wait -n # if there are n_threads_divided iterations running wait here for space to start next job
	fi
done
wait # wait for all jobs in the above loop to be done

echo "Finished aligning ${project_name}" ; date

##################################################################################################
# 7. Launch QC script
##################################################################################################

if [ $run_qc = 'TRUE' ]
then
	echo "launching QC script ${SLURM_SUBMIT_DIR}/bih-alignment/exec/alignment_qc_exec.sh"
	bash ${SLURM_SUBMIT_DIR}/bih-alignment/exec/alignment_qc_exec.sh \
		$project_name \
		$organism
fi

echo "Finished aligning and QC on ${project_name}!" ; date

# change permissions of output folders so all group members can read/write
chmod -R 774 $bam_dir
chmod -R 774 $qc_dir

# # move log
# for x in {a..z}
# do
#        	log_name=/fast/work/groups/ag_sanders/projects/benedict/logs/$(date +%Y%m%d){x}_$(project_name}_variant_calling.txt
#        	if [ ! -f "$log_name" ] ; then
#                	break
#        	fi
# done
# mv /fast/work/groups/ag_sanders/projects/benedict/logs/tmp_aln_log.txt \
#        	$log_name

