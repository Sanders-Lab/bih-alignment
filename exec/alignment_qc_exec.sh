#!/bin/bash
# Run QC script on alignment data, this script should be called at the end of alignment_script.sh
echo 'Running QC script' ; date

##################################################################################################
# 1. Set global options
##################################################################################################
printf '\n ### 1. Setting global options ####\n'

# Needs changing these for each new project:
project_name=$1 # the name of the dir in /fast/groups/ag_sanders/work/data containig the reads (which should contain a dir named fastq/)
echo "project_name set to: $project_name"

n_threads=$(nproc) # number of threads given to slurm job, for best results use a multiple of 4 (e.g. 32, 64, etc.)
n_threads_divided=$(expr $n_threads / 4)
echo "using $n_threads threads"

##################################################################################################
# 2. Activate conda environment
##################################################################################################
# commented out as conda env already actiucated when launching this script from initial alignment script
#printf '\n ### 2. Activating conda environment ####\n'

# dependencies: r samtools picard bwa fastqc multiqc bedtools
# the conda env can be installed with: conda env create -f //fast/groups/ag_sanders/work/tools/conda_envs/alignmentenv_20220505.yml

# test if conda env is present
#condaenvexists=$(conda env list | grep 'alignmentenv' | awk '{print $1}')
#[ -z ${condaenvexists} ] && { echo "ERROR: alignmentenv conda env not found in conda env list" ; exit ; } || echo ''
#[ ${condaenvexists} = "alignmentenv" ] && echo "alignmentenv conda env found in conda env list" || { echo "ERROR: alignmentenv conda env not found i
#n conda env list" ; exit ; }

#conda init bash
#conda activate alignmentenv

##################################################################################################
# 3. Initiation
##################################################################################################
printf '\n ### 3. Initiating script #####\n'

# assign/create directories
bam_dir=//fast/groups/ag_sanders/work/data/${project_name}/bam
qc_dir=//fast/groups/ag_sanders/work/data/${project_name}/qc ; mkdir -m 775 $qc_dir
statsdir=${qc_dir}/alignment_stats ; mkdir -m 775 $statsdir

sortedwc=$(ls $bam_dir | grep 'sorted' | wc -l)
[ $sortedwc = 0 ] && sortorsorted=sort || sortorsorted=sorted


##################################################################################################
# 4. Run mosaicatcher to calc background & plot Strand-Seq overviews
##################################################################################################
printf '\n ### 4. Running mosaicatcher  #####\n'

moscatchdir=${qc_dir}/mosaicatcher ; mkdir -m 775 $moscatchdir

# Run mosaicatcher in singularity
echo 'launching singulairty from docker://smei/mosaicatcher-pipeline-rpe1-chr3'
#//fast/groups/ag_sanders/work/tools/mosaicatcher/build/mosaic count \ # only works if you have mosaic working in your environment (singularity is far easier)
singularity exec --bind /fast docker://smei/mosaicatcher-pipeline-rpe1-chr3 \
	mosaic count \
	-o ${moscatchdir}/counts.txt.gz -i ${moscatchdir}/counts.info \
    -x //fast/work/groups/ag_sanders/data/reference/exclude/GRCh38_full_analysis_set_plus_decoy_hla.exclude \
    -w 200000 $(ls ${bam_dir}/*.bam)

# run R script to generate mosaicatcher plots
echo 'running mosaicather R script'
Rscript ${SLURM_SUBMIT_DIR}/bih-alignment/exec/mosaicatcher_qc.R \
	${moscatchdir}/counts.txt.gz \
	${moscatchdir}/counts.info \
	${moscatchdir}/counts.pdf

echo "mosaicatcher complete, overview Strand-Seq QC plots can be found at: ${moscatchdir}/counts.pdf"

##################################################################################################
# 5. Generate alignment stats
##################################################################################################
printf '\n ### 5. Generate alignement stats #####\n'

flagstats_dir=${statsdir}/flagstats ; mkdir -m 775 $flagstats_dir
idxstats_dir=${statsdir}/idxstats ; mkdir -m 775 $idxstats_dir
bychromdp_dir=${statsdir}/meandepthbychrom ; mkdir -p -m 775 $bychromdp_dir
binneddp_dir=//fast/groups/ag_sanders/scratch/${project_name}/depthstats ; mkdir -p -m 775 $binneddp_dir
insert_dir=${statsdir}/insertsizes ; mkdir -m 775 $insert_dir
insert_hist_dir=${insert_dir}/histograms ; mkdir -m 775 $insert_hist_dir
insert_samps_dir=${insert_dir}/samples ; mkdir -m 775 $insert_samps_dir

libraries=$(ls ${bam_dir}/*bam | sed 's/'.${sortorsorted}.mdup.bam'//g' | sed 's?'${bam_dir}/'??g')

for library in $libraries; do
	(
	echo Generating alignment statistics on ${library}

	bamfile=${bam_dir}/${library}.${sortorsorted}.mdup.bam

	# flagstats
	samtools flagstats -@ 4 $bamfile > ${flagstats_dir}/${library}_flagstats.txt

	# alignment statisitcs per chromosome
	samtools idxstats -@ 4 $bamfile > ${idxstats_dir}/${library}_idxstats.txt

	# depth stats per 200kb bin
	# bedtools coverage -mean \
 #        -a //fast/groups/ag_sanders/work/data/reference/bedfiles/hg38_bedfile_200kb_intervals.bed \
 #        -b $bamfile > ${binneddp_dir}/${library}_mean_depth_200kb_bin.txt

	# depth stats per chromo
	bedtools coverage -mean \
        -a //fast/groups/ag_sanders/work/data/reference/bedfiles/hg38_chromosome_lengths.bed \
        -b $bamfile > ${bychromdp_dir}/${library}_mean_depth_bychrom.txt

    # insert sizes
    picard CollectInsertSizeMetrics -I $bamfile -O ${insert_samps_dir}/${library}_insertsizes.txt \
	   -H ${insert_hist_dir}/${library}_insertsizes.pdf \
	   --QUIET true --VERBOSITY ERROR # makes log easier to read
	) &
	if [[ $(jobs -r -p | wc -l) -ge $n_threads_divided ]]; # allows n_threads number of jobs to be executed in parallel
	then
		wait -n # if there are n_threads jobs running wait here for space to start next job
	fi
done
wait # wait for all jobs in the above loop to be done

gzip -f ${binneddp_dir}/*
gzip -f ${insert_samps_dir}/*

##################################################################################################
# 5. Extract QC stats for plotting
##################################################################################################
printf '\n ### 5. Extract QC stats for plotting  #####\n'

echo 'library gc_content n_reads n_reads_mapped mean_insert_size' | tr " " "\t" > ${statsdir}/all_samples_qc_metrics.txt

for library in $libraries; do
	(
	bamfile=${bam_dir}/${library}.${sortorsorted}.mdup.bam

	# no. of GC bases in BAM file
	gc_ln=$(samtools view -@ 4 $bamfile | awk -F'\t' '{print $10}' | tr -d '\n' | grep -E -o "G|C" | wc -l)
    # total no. of bases in BAM file
	all_ln=$(samtools view -@ 4 $bamfile | awk -F'\t' '{print $10}' | tr "\n" " " | sed 's/\s\+//g' | wc -m)
    # GC content calc
	gc_content=$(echo $gc_ln / $all_ln | bc -l)
	[ -z "$gc_content" ] && gc_content="NA"

    # No. of reads in BAM
	n_reads=$(samtools view -@ 4 $bamfile | wc -l)
	[ -z "$n_reads" ] && n_reads="NA"
	# No. of reads mapped to primary 24 chromosomes
	n_reads_mapped=$(head -n24 ${idxstats_dir}/${library}_idxstats.txt | awk '{print $3}' | paste -sd+ | bc)
	[ -z "$n_reads_mapped" ] && n_reads_mapped="NA"

	# mean insert size
	mean_insert=$(zcat ${insert_samps_dir}/${library}_insertsizes.txt.gz | head -n8 | tail -n1 | awk '{print $6}')
	[ -z "$mean_insert" ] && mean_insert="NA"

    # save output to file
	echo $library $gc_content $n_reads $n_reads_mapped $mean_insert | tr " " "\t" >> ${statsdir}/all_samples_qc_metrics.txt
	) &
	if [[ $(jobs -r -p | wc -l) -ge $n_threads_divided ]]; # allows n_threads number of jobs to be executed in parallel
   	then
       	wait -n
   	fi
done
wait # wait for all jobs in the above loop to be done

##################################################################################################
# 7. running R script to plot QC metrics
##################################################################################################
printf '\n ### 7. running R script to plot QC metrics  #####\n'

Rscript ${SLURM_SUBMIT_DIR}/bih-alignment/exec/alignment_qc_exec.R \
	$project_name \
	$n_threads

echo "Finished alignment script on ${project_name}!" ; date

