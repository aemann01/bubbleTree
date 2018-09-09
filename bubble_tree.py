#!/usr/bin/env python3

'''This script draws a normalized sample abundance heatmap or bubble chart ordered by the leaves on a phylogenetic tree. 
Samples are ordered and colored by a metadata category. Normalized across the row from 0 to 1. 
Useage: python3 bubble_tree.py -i biom.txt -m mapping.txt -t tree.newick -c category -d [heatmap|bubblechart] <-r False -p False>'''

##TO DO: ADD OPTIONS FOR DIFFERENT TREE FORMATS
##TO DO: ALIGN TIP NAMES WITH DOTTED LINES NATIVELY
##TO DO: INLCUDE NATIVE ALIGNMENT/DISTANCE MATRIX BUILD?

import argparse
parser = argparse.ArgumentParser()
requireparser = parser.add_argument_group('required arguments')

requireparser.add_argument('-i', '--input', help='Absolute abundance biom table. Must be tsv formatted.', required=True)
requireparser.add_argument('-t', '--tree', help='Newick formatted tree', required=True)
requireparser.add_argument('-m', '--map', help='Mapping file with metadata corresponding to samples. Must be tsv formatted', required=True)
requireparser.add_argument('-c', '--category', help='Column category from mapping file to order/color samples by', required=True)
requireparser.add_argument('-d', '--display', help='Display data as a heatmap or bubblechart', default="bubblechart")

parser.add_argument('-f', '--treeformat', help='Optional: set tree format type. Default is newick formatted tree.', default='newick')
parser.add_argument('-r', '--remote', help='Set this option as True running on a remote cluster. Disables the automatic $DISPLAY environment varible used by matplotlib', type=bool, default='False')
parser.add_argument('-p', '--previewtree', help='Set this option as True if you want to preview an ASCII version of the imported tree in the standard output', type=bool)

args = parser.parse_args()

if args.remote is not None: 
	import matplotlib
	matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import pandas as pd
import networkx as nx
import numpy as np
import seaborn as sns
import os
from Bio import Phylo
from itertools import repeat

biom = pd.read_csv(args.input, sep="\t", skiprows=1) #load in tab separated biom table -- must have hashed out first row
#print("Reading %s as tsv formatted biom file..." % args.input)
tree = Phylo.read(args.tree, args.treeformat) #load in newick format tree
#print("Reading %s as %s formatted tree file..." % (args.tree, args.treeformat))
metadat = pd.read_csv(args.map, sep="\t") #load in mapping file
#print("Reading in %s as mapping file...\n" % args.map)
#print("Generating figure for %s metadata category..." % args.category)

def normbubble(x,rmin,rmax):
	xnorm = (x - rmin)/(rmax - rmin)
	return xnorm*100

def normheat(x,rmin,rmax):
	xnorm = (x - rmin)/(rmax - rmin)
	return xnorm

def gennormbubble():
	biom['rowmax'] = biom.max(axis=1) #first get max and min value for each row and append to dataframe
	biom['rowmin'] = biom.min(axis=1)
	normdf = pd.DataFrame()
	samps = list(biom['#OTU ID'][1:]) #save first column to append to norm dataframe
	for i in range(1, len(biom)):
		normdata = biom.iloc[i][1:-2].apply(normbubble, args=(biom['rowmin'][i], biom['rowmax'][i])).fillna(0)
		normdf = normdf.append(normdata)
	normdf.insert(loc=0, column="#OTU ID", value=samps)
	reorderbiom(normdf)

def gennormheat():
	biom['rowmax'] = biom.max(axis=1) #first get max and min value for each row and append to dataframe
	biom['rowmin'] = biom.min(axis=1)
	normdf = pd.DataFrame()
	samps = list(biom['#OTU ID'][1:]) #save first column to append to norm dataframe
	for i in range(1, len(biom)):
		normdata = biom.iloc[i][1:-2].apply(normheat, args=(biom['rowmin'][i], biom['rowmax'][i])).fillna(0)
		normdf = normdf.append(normdata)
	normdf.insert(loc=0, column="#OTU ID", value=samps)
	reorderbiom(normdf)

def reorderbiom(normdf):
	leaves = [] #get order of leaves from tree
	for leaf in tree.get_terminals():
		leaves.append(leaf.name)
	subset = normdf.loc[normdf['#OTU ID'].isin(leaves)].set_index('#OTU ID') #pull nodes from biom and reorder by leaves in tree
	ordered = subset.reindex(leaves)
	#group samples by metadata category in mapping file
	filtmeta = metadat[metadat['#SampleID'].isin(list(ordered.columns))] #first remove rows that are not in the biom file
	grouped = filtmeta.groupby(args.category)['#SampleID'].apply(list)
	sampOrder = []
	for i in grouped:
		sampOrder += i
	final = ordered.reindex(columns=sampOrder)
	final.dropna(axis=1, how='all', inplace=True) 	#remove samples if full column is na (i.e., the sample is not present in the biom file)
	if args.display == 'heatmap':
		heat(final)
	elif args.display == 'bubblechart':
		bubble(final, grouped, filtmeta)

def bubble(final, grouped, filtmeta):
	x,y = np.meshgrid(final.columns, final.index) #flatten the matrix
	colmap = sns.color_palette("Set2") + sns.color_palette("Paired") + sns.color_palette("dark") #get color scheme based on grouping category (can do up to 30 categories, add more palettes for more)
	col = []
	j = 0
	for i in grouped:
		col.extend(repeat(colmap[j], len(i)))
		j += 1
	legendCol = sorted(set(col), key=lambda x: col.index(x)) #get corresponding legend values (preserve order in color list)
	legendName = list(filtmeta.groupby(args.category).groups.keys())
	legendGen = [] #generate legend
	j = 0
	for i in legendName:
		legendGen.append(mpatches.Patch(color=legendCol[j], label=legendName[j]))
		j += 1
	#set up subplot aesthetics
	gs = gridspec.GridSpec(1, 2, width_ratios=[0.5, 3]) 
	gs.update(wspace=0.25, hspace=2)
	treeax=plt.subplot(gs[0], frame_on=False)
	bubbleax = plt.subplot(gs[1])
	plt.rc('font', size=0)

	Phylo.draw(tree, axes=treeax, do_show=False)
	plt.scatter(x=x.flatten(), y=y.flatten(), s=final.values.flatten(), zorder=3, c=col,edgecolors="black", axes=bubbleax)


	ax = plt.gca()
	ax.set_ylim(ax.get_ylim()[::-1]) #flip y axis to match tree
	ax.grid(True, linestyle="dotted", linewidth=0.2)
	plt.legend(handles=legendGen, prop={'size': 6}) 
	plt.xticks(rotation=90)
	ax.tick_params(axis = 'both', which = 'major', labelsize = 5)
	treeax.set_ylabel('')
	treeax.set_xlabel('')
	treeax.set_xticks([])
	treeax.set_yticks([])
	#save
	#plt.show() #turn on for testing only
	plt.savefig('%s_bubblePlot.pdf' % args.category, bbox_inches='tight')

##TO DO: MAKE THIS PRETTY SON
def heat(final):
		#set up plot aesthetics
	gs = gridspec.GridSpec(1, 2, width_ratios=[0.5, 3]) 
	gs.update(wspace=0.25, hspace=2)
	treeax=plt.subplot(gs[0], frame_on=False)
	bubbleax = plt.subplot(gs[1])
	ax = plt.gca()
	ax.set_ylim(ax.get_ylim()[::-1]) #flip y axis to match tree
	ax.tick_params(axis = 'both', which = 'major', labelsize = 5)
	bubbleax.set_ylabel('')
	#plot
	plt.rc('font', size=-0)
	treeax.set_ylabel('')
	treeax.set_xlabel('')
	treeax.set_xticks([])
	treeax.set_yticks([])
	Phylo.draw(tree, axes=treeax, do_show=False)
	sns.set(font_scale=1)
	sns.heatmap(final, mask=False, cmap="YlGnBu")
	#save
	#plt.show() #turn on for testing only
	plt.savefig('%s_heatPlot.pdf' % args.category, bbox_inches='tight')

def main():
	assert os.path.exists(args.input), 'Error! File does not exist: %s. Is the path correct?' % args.input
	assert os.path.exists(args.tree), 'Error! File does not exist: %s. Is the path correct?' % args.tree
	assert os.path.exists(args.map), 'Error! File does not exist: %s. Is the path correct?' % args.map
	#preview tree topology?
	if args.previewtree is not None:
		print("Tree preview:\n")
		Phylo.draw_ascii(tree)
	#heatmap or bubble plot?
	if args.display == "heatmap":
		gennormheat()
	elif args.display == "bubblechart":
		gennormbubble()
	else: 
		print("Need to set display parameter! Either heatmap or bubblechart")
main()


