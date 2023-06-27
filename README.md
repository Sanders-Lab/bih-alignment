# bih-alignment

Scripts for performing alignment of paired-end NGS reads to the hg38 reference genome using BWA-MEM on the BIH HPC cluster. 

## Download raw sequencing data

You can download raw fastq data from the genomics facility file exchange server using your MDC account (by changing user@MDC-BERLIN) as follows (adapt link from Email):

```
wget -rnd -np --user=user@MDC-BERLIN --ask-password -A '*.gz'  https://file-exchange.bihealth.org/c2eda02f-523b-4d1a-8909-aea0cd0f5a2d/
```

Find more information on how to get raw data here: https://bimsbstatic.mdc-berlin.de/genomics/howto/user_transfer_fileboxes.html.
Note that file boxes will be kept for 14 days only.

Depending on how closely the Genomics Platform paid attention to our submission form, you may have to edit the filenames of the FASTQ files to get them into our desired format (i.e. `P3069_i301.1.fastq.gz`). 
This can be achieved with the following script: `/fast/groups/ag_sanders/work/projects/benedict/fastq_name_convert/fastq_name_convert.sh`, which compares the indices in the FASTQ header to those that we submitted , and uses this information to rename the file accordingly. 

## Installation

You can download this repository like so:

```
git clone https://github.com/Sanders-Lab/bih-alignment
```

And install the required conda environment (which has all required packages) into your BIH cluster conda workspace like so (you only need to do this once and can then re-use it):

```
conda env create --force --file bih-alignment/alignmentenv_20220905.yml
```

## Usage 

Once the repo is cloned you can launch the complete alignment and QC pipeline like so:
```
sbatch \
    -J alignment \
    -o /fast/work/groups/ag_sanders/projects/${myname}/logs/$(date +%Y%m%d)_${project_name}_alignment.txt \
    bih-alignment/scripts/alignment_script.sh \
    $project_name \
    .1.fastq.gz \
    human
```
* Where `$project_name` is the the name of the directory in `/fast/groups/ag_sanders/work/data` containing the reads, which should contain a dir named `fastq/` with the read files (e.g. set to `P1593` to align reads in `/fast/groups/ag_sanders/work/data/P1593/fastq`). 

* The second command line variable (`.1.fastq.gz` in this example) should be the shared suffix of the first mate FASTQ files. The suffix of the second mate is assumed to be the same with a 2 in place of the 1.

* The third command line variable is the organism to which you want to align the reference genome, currently this can either be set to `human` or `mouse`.

You can either edit `${myname}` manually, or (to make your life easier) add them as an environmental bash variables (e.g. by adding `export myname=benedict` to your `~/.bashrc` file). 

To ensure the pipeline can locate and execute the scripts in the `exec/` directory, please submit the slurm job from the directory into which you initially cloned this repository (i.e. run `sbatch bih-alignment/scripts/alignment_script.sh` as in the example above). 

This repo contains the following scripts:

* `alignment_script.sh` is the complete pipeline and is therefore the recommended script to use.

* `alignment_qc_script.sh` can be run if you only wish to run quality control on data that is already aligned.

The files in `exec/` are called by the main scripts and should not be executed in isolation.  

## Usage (QC only)

If you already have aligned BAM files and wish to run the quality control script in isolation, you can do so like this:
```
sbatch \
    -J alnQC \
    -o /fast/work/groups/ag_sanders/projects/${myname}/logs/$(date +%Y%m%d)_${project_name}_alignment_qc.txt \
    bih-alignment/scripts/alignment_qc.sh \
    $project_name \
    human
```


## Configuration (Optional)

The following parameters can be changed in the global options section in `scripts/alignment_scripts.sh` at the start of each new project:

* `project_name`: the name of the folder in `//fast/groups/ag_sanders/work/data` containing the reads (which should contain a dir named fastq/). Defaults to command line argument 1.

* `memperthread` the memory assigned per thread in slurm job (recommended 2-4G)

* `run_qc`: whether the alignment QC script should be run automatically after alignment is complete [TRUE/FALSE] (default = TRUE)

With these parameters set correctly you should be able to run the pipeline. Further information on each step of the pipeline can be found in the comments of each script.

## Authors 

Code written by Benedict Monteiro & Patrick Weidner. Please contact us with any problems or submit them as an issue in this Github repository.
