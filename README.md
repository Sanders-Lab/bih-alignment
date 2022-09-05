# bih-alignment

Scripts for performing alignment of paired-end NGS reads to the hg38 reference genome using BWA-MEM on the BIH HPC cluster. 

## Installation

You can download this repository like so:

```
git clone https://github.com/benedict909/bih-alignment
```

And install the required conda environment (which has all required packages) into your BIH cluster conda workspace like so:

```
conda env create --force --file bih-alignment/alignmentenv_20220905.yml
```

## Usage 

Once the repo is cloned you can launch the complete alignment and QC pipeline like so:
```
sbatch \
  -J alignment \
  -o /fast/work/groups/ag_sanders/projects/${myname}/logs/$(date +%Y%m%d)_${project_name}_alignment.txt \
  --mail-user=${myemail} \
  bih-alignment/scripts/alignment_script.sh \
  $project_name
```
Where `$project_name` is the the name of the directory in `/fast/groups/ag_sanders/work/data` containing the reads, which should contain a dir named `fastq/` with the read files (e.g. set to `P1593` to align reads in `/fast/groups/ag_sanders/work/data/P1593/fastq`). 

You can either edit `${myemail}` and `${myname}` manually, or (to make your life easier) add them as an environmental bash variables (e.g. by adding `export myname=benedict` to your `~/.bashrc` file). 

This repo contains the following scripts:

* `alignment_script.sh` is the complete pipeline and is therefore the recommended script to use.

* `alignment_qc_script.sh` can be run if you only wish to run quality control on data that is already aligned.

The files in `/exec` are called by the main scripts and should not be executed in isolation.  


## Configuration (Optional)

The following parameters can be changed in the global options section in `scripts/alignment_scripts.sh` at the start of each new project:

* `project_name`: the name of the folder in `//fast/groups/ag_sanders/work/data` containing the reads (which should contain a dir named fastq/). Defaults to command line argument 1.

* `memperthread` the memory assigned per thread in slurm job (recommended 2-4G)

* `mate1_suffix`: the suffix of read pair mate 1 fastq file (e.g. `_R1.fastq.gz`)

* `run_qc`: whether the alignment QC script should be run automatically after alignment is complete [TRUE/FALSE] (default = TRUE)

With these 4 parameters set correctly you should be able to run the pipeline. Further information on each step of the pipeline can be found in the comments of each script.
