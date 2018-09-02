import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from Bio import Phylo

#load in tab separated biom table -- must have hashed out first row
biom = pd.read_csv("test_biom.txt", sep="\t", skiprows=1)

#load in newick format tree
tree = Phylo.read("test_newick.tre", "newick")

#preview tree
print("Tree preview:\n")
Phylo.draw_ascii(tree)

#get order of leaves from tree
leaves = []
for leaf in tree.get_terminals():
	leaves.append(leaf.name)

#pull nodes from biom and reorder by leaves in tree
subset = biom.loc[biom['#OTU ID'].isin(leaves)].set_index('#OTU ID')
ordered = subset.reindex(leaves)
#flatten the matrix to plot cleveland style dotplot
x,y = np.meshgrid(ordered.columns, ordered.index)

#plot the bubble chart, save to pdf
#plt.subplot(1,2,2)
plt.scatter(x=x.flatten(), y=y.flatten(), s=ordered.values.flatten())
#flip y axis to match tree
ax = plt.gca()
ax.set_ylim(ax.get_ylim()[::-1])
plt.savefig("bubble.pdf")

#plot the tree, save to pdf
#plt.subplot(1,2,1)
Phylo.draw(tree)
plt.savefig("tree.pdf")
