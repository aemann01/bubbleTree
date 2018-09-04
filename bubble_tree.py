import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import seaborn as sns
import matplotlib.patches as mpatches
from Bio import Phylo
from itertools import repeat

##ADD ARGUMENT PARSER HERE

#load in tab separated biom table -- must have hashed out first row
biom = pd.read_csv("test_biom.txt", sep="\t", skiprows=1)

#load in newick format tree
tree = Phylo.read("test_newick.tre", "newick")

#load in mapping file
metadat = pd.read_csv("test_map.txt", sep="\t")

#preview tree
print("Tree preview:\n")
##print("Number of leaves %i" % //////)

Phylo.draw_ascii(tree)

#get order of leaves from tree
leaves = []
for leaf in tree.get_terminals():
	leaves.append(leaf.name)

#pull nodes from biom and reorder by leaves in tree
subset = biom.loc[biom['#OTU ID'].isin(leaves)].set_index('#OTU ID')
ordered = subset.reindex(leaves)

#####EXTRA GROUPING/COLOR OPTIONS
#group samples by metadata category in mapping file
grouped = metadat.groupby('CatC')['SampleID'].apply(list)
sampOrder = []
for i in grouped:
	sampOrder += i
final = ordered.reindex(columns=sampOrder)

#flatten the matrix to plot cleveland style dotplot
x,y = np.meshgrid(final.columns, final.index)

#plot the bubble chart, save to pdf
#plt.subplot(1,2,2)
#get color scheme based on grouping category
colmap = sns.color_palette()

col = []
j = 0
for i in grouped:
	col.extend(repeat(colmap[j], len(i)))
	j += 1

#get corresponding legend values (preserve order in color list)
legendCol = sorted(set(col), key=lambda x: col.index(x)) 
legendName = metadat.CatC.unique()

#generate legend
legendGen = []
j = 0
for i in legendCol:
	legendGen.append(mpatches.Patch(color=legendCol[j], label=legendName[j]))
	j += 1

plt.scatter(x=x.flatten(), y=y.flatten(), s=final.values.flatten(), zorder=3, c=col,edgecolors="black")
#flip y axis to match tree
ax = plt.gca()
ax.set_ylim(ax.get_ylim()[::-1])
ax.grid(True)
plt.legend(handles=legendGen)
plt.xticks(rotation=90)
plt.savefig("bubble.pdf")

#plot the tree, save to pdf
#plt.subplot(1,2,1)
Phylo.draw(tree)
plt.savefig("tree.pdf")








