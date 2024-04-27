# SSN

## Files
- `ssn.zip` - zip archive containing the input network as a file of edges
- `cytoscape_sessiions.zip` - zip archive containing cytoscape sessions of the annotated SSN
    - `ssn_annotated_characterized_types.cys` - colored characterized TPSs according to their type
    - `ssn_annotated_superkingdoms.cys` - colored sequences according to the superkingdoms
    - `ssn_annotated_novelty.cys` - colored sequences according to the novelty score
    - `ssn_annotated_reliability.cys` - colored sequences according to the reliability score

## Cytoscape
Cytoscape can be downloaded here: https://cytoscape.org/download.html

The Cytoscape version used in this project was: 3.10.1
```
Version: 3.10.1 
Java: 17.0.5 by Eclipse Adoptium
OS: Windows 10 10.0 - amd64
```

To use Cytocape, see the Cytoscape manual: https://manual.cytoscape.org/en/stable/

### py4cytoscape 
`py4cystocape` is a Python package to communicate with Cytoscape, see: https://py4cytoscape.readthedocs.io/en/latest/

See an overview jupyter notebook here: https://github.com/cytoscape/cytoscape-automation/blob/master/for-scripters/Python/Overview-of-py4cytoscape.ipynb

Useful case for this project was for example to programatically export images of coloring by different domains, see `../notebooks/SSN/01_color_domains.ipynb`