# Phylogenetic tree

## Files
- `tps_mining.tree` - phylogenetic tree in Newick format constructed with Fasttree2 from representatives of 50% sequence identity clusters from the mining and from the characterized TPSs from TPS db
- `itol_annotations` - this directory contains annotations which were loaded in iTOL to annotate the tree

## iTOL
iTOL (Interactive Tree of Life) is an online tool to display and annotate phylogenetic trees, see: https://itol.embl.de/

See iTOL manual page here: https://itol.embl.de/help.cgi

For this project, following steps were done:

The tree can be found here: https://itol.embl.de/tree/91113129162299121709309485 (or upload the tree through Upload menu)


Set the view through the control panel. 

I used:
- Basic:
    - Mode: Circular
    - Branch lengths: Ignore
    - Labels: Hide
- Advanced:
    - Node options: Collapesed nodes -> Un-collapse all

Upload the annotations by drag and drop of the annotation files in the `itol_annotations` directory.

With a free iTOL account, the tree is saved without the annotations but they can be always easily loaded by this drag and drop of annotation files. 

![SSN](https://github.com/CalounovaT[TPS_mining]/blob/main/03_mining_analysis/plots/annotated_tree.png?raw=true)
