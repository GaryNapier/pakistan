import sys
import argparse
from Bio import Phylo
from ete3 import Tree
import re
import os

def main(args):

    infile = args.infile
    outfile = args.outfile
    year = int(args.cut)

    tree_file = args.infile
    # tree_file_newick = tree_file + ".nwk"

    # Convert and write Beast MCC nexus file to newick so ete can read
    # Phylo.convert(tree_file, "nexus", tree_file_newick, "newick")

    # Read in newick tree and parse as tree with ete
    # ETE - see http://etetoolkit.org/docs/latest/tutorial/tutorial_trees.html
    # t = Tree(tree_file_newick, format = 1)
    t = Tree(tree_file, format = 1)

    # Get all sample names
    samps = t.get_leaf_names()

    # Get years of samples
    years = []
    for samp in samps:
        samp_split = samp.split("_")
        years.append(samp_split[-1])
    # Get max date
    max_date = float(max(years))

    # Add node names if none:
    for i, node in enumerate(t.traverse("preorder")):
        if node.name == '':
            node.name = node.name + str(i+1)

    # Get all dists from root to get max dist
    dists_from_root = []
    for i, node in enumerate(t.traverse("preorder")):
        if node.is_root():
            rt = t&node.name
        dist = rt.get_distance(node.name)
        dists_from_root.append(dist)
    max_dist = max(dists_from_root)

    # Loop over tree and convert dists from root to year, adding correction
    # Add the sample names (leaves) to the node_list if node they belong to is at the cutoff point and the node is not a single sample
    node_list = []
    for i, node in enumerate(t.traverse("preorder")):
        # Convert to year - have to add this correction (dist from root - max dist of all distances in tree) + max date
        # See ptreeFromPhylo function in TransPhylo package - https://github.com/xavierdidelot/TransPhylo/blob/master/R/ptreeFromPhylo.R
        year_dist = (rt.get_distance(node.name) - max_dist) + max_date
        if year_dist >= max_date - year and len(node) > 1:
            leaves = node.get_leaf_names()
            node_list.append(leaves)

    # Remove clades that are subsets of higher clades
    final_node_list = node_list[:]
    for node in node_list:
        for node_again in node_list:
            if set(node).issubset(set(node_again)) and node != node_again:
                final_node_list.remove(node)
                break

    # Set up table to store clusters (to import to R)
    clust_dict = {}
    for i, node in enumerate(final_node_list):
        clust_dict[str(i+1)] = node

    # Write out the dict as table
    with open(outfile,"w") as O:
        for key in clust_dict:
            for val in clust_dict[key]:
                O.write(",".join([key,val]) + "\n")
    O.close()


parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--infile',type=str, default='', help='Name of BEAST MCC tree file', required=True)
parser.add_argument('--outfile',type=str, default='', help='Name of clusters output file', required=True)
parser.add_argument('--cut',type=str, default='', help='Number of years before last date at which to cut tree. e.g. if date of last sample is 2000, and want to cut tree at 1990, enter 10')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
