# Computing irreproducible discovery rate (IDR)

Nextflow workflow to calculate IDR values and evaluate self-consistency and consistency between replicates or pseudo-replicates of next-generation sequencing libraries that capture regions of enriched signal, such as ChIP-seq or ATAC-seq.

## Author

Yuvia Alhelí PÉREZ-RICO

Affiliation: European Molecular Biology Laboratory | EMBL

## Dependencies

To use this workflow, you will need to install the following programs, the indicated version is the one that I have used:

- nextflow (20.04.1)
- samtools (1.9)
- bedtools (v2.29.2)
- macs2 (2.2.7.1)
- IDR (2.0.4.2)

Note: I suggest to install all programs in a conda environment.

## How to use the workflow

Thank you for your interest in using the workflow!

After downloading the main script and the configuration file, indicate in the configuration file the type of comparison that you would like to perform (self-consistency, replicates or pseudo-replicates), as well as, values for ranking and peak file format. Then, prepare the following file and change the paths in the configuration file accordingly:

- A comma-separated file indicating sample names and paths to bam or peak files (examples for the three different types of comparisons can be found in the 'docs' folder).

If you are not using SLURM, then do additional modifications to the configuration file considering the workload manager that you use. Finally, write a simple bash script that will be submitted to the cluster to activate the conda environment and start the main nextflow job:

`source /home/user/miniconda2/bin/activate /home/user/conda-envs/IDR`

`nextflow run IDR.nf`

## Acknowledgements

This repository is part of a project that has received funding from the European Union’s Horizon 2020 research and innovation programme under the Marie Skłodowska-Curie grant agreement No 882771.

## License

Licenced under the GNU general public license.

