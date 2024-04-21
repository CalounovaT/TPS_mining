# TPS db analysis

## 1. Create and process files for analysis

Snakemake pipeline to prepare various files for further analysis
- TPS db preprocessing 
- motif scanning
- Pfam and SUPERFAMILY scanning
- Pfam and SUPERFAMILY domain architectures
- prenyltrasnferases similarity
- ESM2 embeddings

See `Snakefile`

## 2. Notebooks

In the `notebooks` directory, find various notebooks analyzing the TPS db

- `00_prepare_df.ipynb`
    - small modifcations to the dataframe
- `01_motif_presence_analysis.ipynb`
    - analysis of the conserved motifs presence
- `02_domain_architecture_analysis.ipynb`
    - analysis of the domain archotectures
- `03_length_analysis.Rmd`
    - analysis of the length distribution
- `04_product_occurences_counts.ipynb`
    - count occurences of all observed products
    - assign a most "rare" product count to each sequence
- `05_embeddings_visualization.ipynb`
    - t-SNE projection of the ESM2 protein embeddings

## 3. Plots

Resulting plots can be found in the `plots` directory