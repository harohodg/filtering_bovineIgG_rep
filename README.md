### Nextflow workflows for characterizing bovine immunological nanopore sequencing data
This workflow has been tested on [Narval](https://docs.alliancecan.ca/wiki/Narval/en) but should work on any of the other [Digital Research Alliance of Canada](https://alliancecan.ca/en) systems. With a bit of editing these scripts should be able to run on any system with [Nextflow](https://www.nextflow.io/docs/latest/index.html) and [Fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [Fastp](https://github.com/OpenGene/fastp) and [Seqkit](https://bioinf.shenwei.m) and [Bioawk](https://github.com/lh3/bioawk) installed. 


Note : The nextflow scripts expects compressed fastq files.
Note: Nextflow doesn't like relative paths for input files. `~/somewhere/or/other` works as does `$(realpath a/relative/path)`

#### Single data-set
In an interative job

```
module load StdEnv/2020 nextflow/23.04.3

nextflow run <path/to/fastqc_and_fastp.nf> --input_file <input_file.fastq.gz> --output_dir <path/to/results_folder>
nextflow run <path/to/IgG_filtering.nf> --input_file <fastp-trimmed-file.fastq.gz> --output_dir <path/to/results_folder>
```



#### Multiple data-sets

1. Check quality and trim sequences using Fastqc and Fastp
Assuming you've put the basecalled files in `~/scratch/bovine_nanopore_data` and that you are currently in the folder that you want all the results in.
```
export SCRIPTS_DIR="path/to/scripts/folder"

module load meta-farm/1.0.2

farm_init.run fastqc_fastp-farm

find ~/scratch/bovine_nanopore_data -name '*.fastq.gz' | parallel --dry-run 'NXF_WORK=$SLURM_TMPDIR/work nextflow run '${SCRIPTS_DIR}'/fastqc_and_fastp.nf --input_file {}  --output_dir '$(pwd)/'$(echo "{/.}" | sed "s/.fastq//")_results' > fastqc_fastp-farm/table.dat

eval cp ${SCRIPTS_DIR}/fastqc_and_fastp-job_script.sh fastqc_fastp-farm/job_script.sh
eval cp ${SCRIPTS_DIR}/single_case.sh fastqc_fastp-farm/single_case.sh

cd fastqc_fastp-farm && submit.run 2
```



2. Run filtering pipeline on each fastp trimmed file.
Assuming you are still in the fastqc_fastp-farm folder, that the meta-farm module is still loaded, and that the previous jobs are finished
```
cd ..
farm_init.run IgG_filtering-farm

find bc*/fastp -name 'bc*trimmed.fastq.gz' | parallel --dry-run 'NXF_WORK=$SLURM_TMPDIR/work nextflow run '${SCRIPTS_DIR}'/IgG_filtering.nf --input_file '$(pwd)'/{}  --output_dir '$(pwd)/'$(echo "{/.}" | sed "s/-trimmed.fastq//")_results' > IgG_filtering-farm/table.dat

eval cp ${SCRIPTS_DIR}/IgG_filtering-job_script.sh IgG_filtering-farm/job_script.sh
eval cp ${SCRIPTS_DIR}/single_case.sh IgG_filtering-farm/single_case.sh

cd IgG_filtering-farm && submit.run 2 
```


Use query.run inside the farm folders to check on the overall status over the jobs
If you have a choice of which def-account_name to use when you submit jobs you can over ride the default one for the meta-farm jobs by submitting with `submit.run 4 '--account def-other_account'`
