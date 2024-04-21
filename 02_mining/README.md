# TPS mining

The mining was performed using a custom Snakemake pipeline, see `Snakefile` using HMMER3, and also includes 5 filtering steps.
- create TPS Pfam db
- create TPS SUPERFAMILY db
- split the databases
- Pfam mining
- SUPERFAMILY mining
- combining the results
- length filtering (1)
- stronger hit in other family (2)
- motif presence (3)
- wrong architecture (4)
- similar to PTs (5)

The main resulting file is: `results/filtering/all_filtered_5_unique_no_stop.fasta`