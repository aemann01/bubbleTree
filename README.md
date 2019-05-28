# Bubble Tree

![example](img/figexample.png)

NOTE: Under construction!!

Bubble tree takes a tab separated frequency table (e.g., a biom table with taxa as rows and samples as columns), a newick, nexus, nexml, phyloxml, or cdao formatted phylogenetic tree, and a tab separated file with sample metadata to create either a Cleveland dot plot--style or heatmap and tree.

### Prerequisites

Bubble tree is written with python 3+ and relies on the following packages:

* [BioPython](https://biopython.org/) 
* [Matplotlib](https://matplotlib.org/)
* [Pandas](https://pandas.pydata.org/)
* [Numpy](http://www.numpy.org/)
* [Seaborn](https://seaborn.pydata.org/)

### Examples
Generate a heatmap figure ordered by the mapping.txt column Habitat
```
bubble_tree.py -i exampes/biom.txt -t examples/tree.tre -m examples/map.txt -c Habitat -d heatmap
```

Generate a bubble chart figure colored and ordered by the mapping.txt column Species and view the tree structure in ASCII format to screen
```
bubble_tree.py -i examples/biom.txt -t examples/tree.tre -m examples/map.txt -c Species -d bubble -p True
```

Generate a bubble chart figure colored and ordered by the mapping.txt column Genus on a remote cluster (Disables the automatic $DISPLAY environment varible used by matplotlib)
```
bubble_tree.py -i examples/biom.txt -t examples/tree.tre -m examples/map.txt -c Genus -d bubble -r True
```

Help and parameter description
```
bubble_tree.py -h
```
