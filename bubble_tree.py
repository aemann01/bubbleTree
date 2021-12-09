#!/usr/bin/env python3

'''This script draws a normalized sample abundance heatmap or bubble chart ordered by the leaves on a phylogenetic tree. 
Samples are ordered and colored by a metadata category. Normalized across the row (single tree tip) from 0 to 1. 
Useage: python3 bubble_tree.py -i biom.txt -m mapping.txt -t tree.newick -c category -a asvid -s sampleid -d <bubblechart> <-r False> <-p False> <-n row>'''

import argparse
parser = argparse.ArgumentParser()
requireparser = parser.add_argument_group('required arguments')
requireparser.add_argument('-i', '--input', help='Frequency table. Must be tsv formatted.', required=True)
requireparser.add_argument('-t', '--tree', help='Phylogenetic tree', required=True)
requireparser.add_argument('-m', '--map', help='Mapping file with metadata corresponding to samples. Must be tsv formatted', required=True)
requireparser.add_argument('-c', '--category', help='Column category from mapping file to order/color samples by', required=True)
requireparser.add_argument('-a', '--asvids', help='Column name in biom table that contains ASV (or OTU or whatever) ids', required=True)
requireparser.add_argument('-s', '--sampleids', help='Column name in mapping file that contains sample IDs', required=True)
parser.add_argument('-d', '--display', help='Display data as a heatmap or bubblechart', default="bubblechart")
parser.add_argument('-f', '--treeformat', help='Optional: set tree format type. Default is newick formatted tree.', default='newick')
parser.add_argument('-r', '--remote', help='Set this option as True running on a remote cluster. Disables the automatic $DISPLAY environment varible used by matplotlib', type=bool, default='False')
parser.add_argument('-p', '--previewtree', help='Set this option as True if you want to preview an ASCII version of the imported tree', type=bool)
parser.add_argument('-n', '--norm', help='Normalize read counts by row (0-1) or log transformation', default='row')
##TODO: Root tree function
#parser.add_argument('-o', '--outgroup', help='Root tree either at midpoint or with named outgroup', default='Null')
args = parser.parse_args()

if args.remote is not None: 
	import matplotlib
	matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import pandas as pd
import numpy as np
import seaborn as sns
import os
from Bio import Phylo
from itertools import repeat

biom = pd.read_csv(args.input, sep="\t") #load in tab separated biom table 
print("Reading %s as tsv formatted biom file..." % args.input)
tree = Phylo.read(args.tree, args.treeformat) #load in tree
print("Reading %s as %s formatted tree file..." % (args.tree, args.treeformat))
metadat = pd.read_csv(args.map, sep="\t") #load in mapping file
print("Reading in %s as mapping file...\n" % args.map)
print("Generating figure for %s metadata category..." % args.category)

def rownormbubble(x,rmin,rmax):
	xnorm = (x - rmin)/(rmax - rmin)
	return xnorm*100

def rownormheat(x,rmin,rmax):
	xnorm = (x - rmin)/(rmax - rmin)
	return xnorm

def gennormbubble():
	biom['rowmax'] = biom.iloc[:, 1:len(biom.columns)].max(axis=1) #first get max and min value for each row and append to dataframe
	biom['rowmin'] = biom.iloc[:, 1:len(biom.columns)].min(axis=1)
	print(biom)
	normdf = pd.DataFrame()
	samps = list(biom[args.asvids][1:]) #save first column to append to norm dataframe
	if args.norm == "row":
		for i in range(1, len(biom)):
			normdata = biom.iloc[i][1:-2].apply(rownormbubble, args=(biom['rowmin'][i], biom['rowmax'][i])).fillna(0)
			normdf = normdf.append(normdata)
		normdf.insert(loc=0, column=args.asvids, value=samps)	
	elif args.norm == "log":
		for i in range(1, len(biom)):
			normdata = np.log10(biom.iloc[i][1:-2].astype(np.float64))
			normdf = normdf.append(normdata)
		normdf.insert(loc=0, column=args.asvids, value=samps)
	reorderbiom(normdf)

def gennormheat():
	biom['rowmax'] = biom.iloc[:, 1:len(biom.columns)].max(axis=1) #first get max and min value for each row and append to dataframe
	biom['rowmin'] = biom.iloc[:, 1:len(biom.columns)].min(axis=1)
	normdf = pd.DataFrame()
	samps = list(biom[args.asvids][1:]) #save first column to append to norm dataframe
	if args.norm == "row":
		for i in range(1, len(biom)):
			normdata = biom.iloc[i][1:-2].apply(rownormbubble, args=(biom['rowmin'][i], biom['rowmax'][i])).fillna(0)
			normdf = normdf.append(normdata)
		normdf.insert(loc=0, column=args.asvids, value=samps)	
	elif args.norm == "log":
		for i in range(1, len(biom)):
			normdata = np.log10(biom.iloc[i][1:-2].astype(np.float64))
			normdf = normdf.append(normdata)
		normdf.insert(loc=0, column=args.asvids, value=samps)
	reorderbiom(normdf)

def reorderbiom(normdf):
	leaves = [] #get order of leaves from tree
	for leaf in tree.get_terminals():
		leaves.append(leaf.name)
	subset = normdf.loc[normdf[args.asvids].isin(leaves)].set_index(args.asvids) #pull nodes from biom and reorder by leaves in tree
	ordered = subset.reindex(leaves)
	filtmeta = metadat[metadat[args.sampleids].isin(list(ordered.columns))] #first remove rows that are not in the biom file
	grouped = filtmeta.groupby(args.category)[args.sampleids].apply(list) #group samples by metadata category
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
	col = col * len(x)
	legendCol = sorted(set(col), key=lambda x: col.index(x)) #get corresponding legend values (preserve order in color list)
	legendName = list(filtmeta.groupby(args.category).groups.keys())
	legendGen = [] #generate legend
	j = 0
	for i in legendName:
		legendGen.append(mpatches.Patch(color=legendCol[j], label=legendName[j]))
		j += 1	
	gs = gridspec.GridSpec(1, 2, width_ratios=[0.5, 3]) #set up subplot aesthetics
	gs.update(wspace=0.25, hspace=2)
	treeax=plt.subplot(gs[0], frame_on=False)
	bubbleax = plt.subplot(gs[1])
	plt.rc('font', size=0)
	Phylo.draw(tree, axes=treeax, do_show=False)
	plt.scatter(x=x.flatten(), y=y.flatten(), s=final.values.flatten(), c=col, zorder=3, edgecolors="black", axes=bubbleax)
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
	#plt.show() #for debugging
	plt.savefig('%s_bubblePlot.pdf' % args.category, bbox_inches='tight') #save

def heat(final):
	gs = gridspec.GridSpec(1, 2, width_ratios=[0.5, 3]) #set up plot aesthetics
	gs.update(wspace=0.25, hspace=2)
	treeax=plt.subplot(gs[0], frame_on=False)
	bubbleax = plt.subplot(gs[1])
	ax = plt.gca()
	ax.set_ylim(ax.get_ylim()[::-1]) #flip y axis to match tree
	ax.tick_params(axis = 'both', which = 'major', labelsize = 5)
	bubbleax.set_ylabel('')
	plt.rc('font', size=-0)
	treeax.set_ylabel('')
	treeax.set_xlabel('')
	treeax.set_xticks([])
	treeax.set_yticks([])
	Phylo.draw(tree, axes=treeax, do_show=False)
	sns.set(font_scale=1)
	sns.heatmap(final, mask=False, cmap="YlGnBu", xticklabels=True, yticklabels=True)
	#plt.show() #for debugging
	plt.savefig('%s_heatPlot.pdf' % args.category, bbox_inches='tight') #save

def main():
	assert os.path.exists(args.input), 'Error! File does not exist: %s. Is the path correct?' % args.input
	assert os.path.exists(args.tree), 'Error! File does not exist: %s. Is the path correct?' % args.tree
	assert os.path.exists(args.map), 'Error! File does not exist: %s. Is the path correct?' % args.map
	if args.previewtree is not None: #preview tree topology?
		print("Raw tree preview:\n")
		Phylo.draw_ascii(tree)
	if args.display == "heatmap": #heatmap or bubble plot?
		gennormheat()
	elif args.display == "bubblechart":
		gennormbubble()

main()

__author__ = "Allison E. Mann"
__license__ = "GPL"
__version__ = "1.0.1"
__email__="allison.e.mann@gmail.com"

