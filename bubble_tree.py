#!/usr/bin/env python3

'''This script draws a sample abundance bubble chart ordered by the leaves on a phylogenetic tree. Samples are ordered and colored by a metadata category. Useage: python3 bubble_tree.py -b <biom.txt> -m <mapping.txt> -t <tree.newick> -c <category>'''

import argparse
parser = argparse.ArgumentParser()
requireparser = parser.add_argument_group('required arguments')

requireparser.add_argument('-b', '--biom', help='Absolute abundance biom table. Must be tsv formatted.', required=True)
requireparser.add_argument('-t', '--tree', help='Newick formatted tree', required=True)
requireparser.add_argument('-m', '--map', help='Mapping file with metadata corresponding to samples. Must be tsv formatted', required=True)
requireparser.add_argument('-c', '--category', help='Column category from mapping file to order/color samples by', required=True)
parser.add_argument('-f', '--treeformat', help='Optional: set tree format type. Default is newick formatted tree.', default='newick')
parser.add_argument('-r', '--remote', help='Set this option as True running on a remote cluster. Disables the automatic $DISPLAY environment varible used by matplotlib', type=bool, default='False')
parser.add_argument('-p', '--previewtree', help='Set this option as True if you want to preview an ASCII version of the imported tree in the standard output', type=bool)

args = parser.parse_args()

if args.remote is not None:
	import matplotlib
	matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pandas as pd
import networkx as nx
import numpy as np
import seaborn as sns
import matplotlib.patches as mpatches
from Bio import Phylo
from itertools import repeat

##TO DO: AT SOME POINT, ADD ALIGNMENT AND DISTANCE MATRIX OPTIONS?
##TO DO: SET DISTANCE BETWEEN TREE AND BUBBLE PLOT SO NO OVERLAP

###LOAD FILES###
#load in tab separated biom table -- must have hashed out first row
biom = pd.read_csv(args.biom, sep="\t", skiprows=1)
print("Reading %s as tsv formatted biom file..." % args.biom)
#load in newick format tree
tree = Phylo.read(args.tree, args.treeformat)
print("Reading %s as %s formatted tree file..." % (args.tree, args.treeformat))
#load in mapping file
metadat = pd.read_csv(args.map, sep="\t")
print("Reading in %s as mapping file...\n" % args.map)
print("Generating figure for %s metadata category..." % args.category)

###NORMALIZE COUNTS IN BIOM TABLE###
#first get max and min value for each row and append to dataframe
biom['rowmax'] = biom.max(axis=1)
biom['rowmin'] = biom.min(axis=1)

def norm(x,rmin,rmax):
	xnorm = (x - rmin)/(rmax - rmin)
	return xnorm*300
normdf = pd.DataFrame()

#save first column to append to norm dataframe
samps = list(biom['#OTU ID'][1:])

for i in range(1, len(biom)):
	normdata = biom.iloc[i][1:-2].apply(norm, args=(biom['rowmin'][i], biom['rowmax'][i])).fillna(0)
	normdf = normdf.append(normdata)
normdf.insert(loc=0, column="#OTU ID", value=samps)

###REORDER NORMALIZED BIOM TABLE BY BRANCH ORDER###
#preview tree topology
if args.previewtree is not None:
	print("Tree preview:\n")
	Phylo.draw_ascii(tree)

#get order of leaves from tree
leaves = []
for leaf in tree.get_terminals():
	leaves.append(leaf.name)

#pull nodes from biom and reorder by leaves in tree
subset = normdf.loc[normdf['#OTU ID'].isin(leaves)].set_index('#OTU ID')
ordered = subset.reindex(leaves)

###GROUP DATA BY METADATA CATEGORY###
#group samples by metadata category in mapping file
#first remove rows that are not in the biom file
filtmeta = metadat[metadat['#SampleID'].isin(list(ordered.columns))]
grouped = filtmetadat.groupby(args.category)['#SampleID'].apply(list)

sampOrder = []

#sort rows by highest to lowest value 


#if you want to group by the listed metadata categories
for i in grouped:
	sampOrder += i
final = ordered.reindex(columns=sampOrder)
#remove samples if full column is na (i.e., the sample is not present in the biom file)
final.dropna(axis=1, how='all', inplace=True)

###GENERATE CLEVELAND STYLE DOTPLOT AND PLOT###
#flatten the matrix
x,y = np.meshgrid(final.columns, final.index)

#get color scheme based on grouping category
colmap = sns.color_palette("Paired") + sns.color_palette("dark") + sns.color_palette("Set2")

col = []
j = 0
for i in grouped:
	col.extend(repeat(colmap[j], len(i)))
	j += 1

#get corresponding legend values (preserve order in color list)
legendCol = sorted(set(col), key=lambda x: col.index(x)) 
legendName = list(filtmeta.groupby(args.category).groups.keys())

#generate legend
legendGen = []
j = 0
for i in legendName:
	legendGen.append(mpatches.Patch(color=legendCol[j], label=legendName[j]))
	j += 1

gs = gridspec.GridSpec(1, 2, width_ratios=[1, 3]) 
gs.update(wspace=0.5, hspace=2)
treeax=plt.subplot(gs[0], frame_on=False)
bubbleax = plt.subplot(gs[1])

Phylo.draw(tree, axes=treeax, do_show=False)

plt.scatter(x=x.flatten(), y=y.flatten(), s=final.values.flatten(), zorder=3, c=col,edgecolors="black", axes=bubbleax)
ax = plt.gca()
ax.set_ylim(ax.get_ylim()[::-1]) #flip y axis to match tree
ax.grid(True, linestyle="dotted", linewidth=0.2)
plt.legend(handles=legendGen, prop={'size': 6}) 
plt.xticks(rotation=90)
ax.tick_params(axis = 'both', which = 'major', labelsize = 7.5)
treeax.set_ylabel('')
treeax.set_xlabel('')
treeax.set_xticks([])
treeax.set_yticks([])

plt.show() #turn on for testing only
plt.savefig('%s_bubblePlot.pdf' % args.category, bbox_inches='tight')
metadat.close()
biom.close()
tree.close()




