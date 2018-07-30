# Nextflow ASE Pipeline

A Nextflow pipeline to use samtools and WASP to get ASE summary data.

## Usage:

1) Set the variables within the ase_pipeline.nf file in the Configuration for the specific run section, and anything needed for the pipeline.

2) External programs are not included, but can be downloaded publicly and added to a "/bin" folder.  They should be:
-WASP
-your aligner
-Samtools

3) Edit the ase_pipeline.nf process 'remap' to match your aligner and the arguments for the aligner.
