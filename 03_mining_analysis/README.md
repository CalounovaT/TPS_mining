# Mining analysis

## Create and process files for analysis
Snakemake pipeline to prepare various files for further analysis
- sequence similarity of mined sequences to characterized sequences
- 50%SI clustering
- SSN using blast
- phylogenetic tree construction
- ESM2 embeddings
- (SSN using blast for plants only)

## Sequence annotation
Combination of manual and computational annotation for mining results from individual databases is documented in `annotation/{database}/{database}_mapping.ipynb`

## Notebooks
In the `notebooks` directory, find various notebooks analyzing the TPS db
- TODO: describe the individual notebooks
- `00_merge_annotations.ipynb`
    - merge individually annotated mined database files
- `01_domain_analysis.ipynb`
    - basic analysis - Pfam/SUPERFAMILY domain hits and architectures
- `02_taxonomy_analysis.ipynb`
    - analysis of taxonomy
- `03_length_analysis.Rmd`
    - analysis of length distribution

### Cluster analysis
- `01_cluster_analysis.ipynb`
    - basic analysis of clusters from 50%SI clustering (cluster sizes, average identity, superkingdom composition, domain composition, domain architectures, source databases composition)
### Embeddings
- `01_esm2_projections.ipynb`
    - PCA and UMAP projections of the ESM2 embeddings
### Novelty score
- `01_SSN_community_detection.ipynb`
    - automatic cluster detection using Louvain method in networkx, visualization of the clusters in Cytoscape (py4cytoscape)
- `02_taxonomic_novelty.ipynb`
    - creation of taxonomic score based on the lineage of the sample origin and taxomic representation in TPSdb
- `03_novelty_score.ipynb`
    - creation of the novelty score and testing on test set from TPSdb
### Phylogenetics
- `01_analyze_min_distances.ipynb`
    - analyze distances to closest characterized TPSs
- `02_unchar_clade_size.ipynb`
    - get size of clade consisting of only uncharacterized sequences for each mined sequence
- `03_create_annotations.ipynb`
    - create dataset annotations for visualization in iTOL
### Reliability score
- `01_basic_reliability_score.ipynb`
    - definition of reliability score
### SSN
- `01_color_domains.ipynb`
    - automaticaly color each domain in Cytoscape and generate pictures
- `02_product_uniqness.ipynb`
    - exploration of rarest product count in cytoscape

## Plots
Resulting plots can be found in the `plots` directory
