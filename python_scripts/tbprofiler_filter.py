#! /usr/bin/env python

# tbprofiler_filter.py --metadata <metadata file> --clusters-file <clusters file> --tbp-results <tbprofiler results directory> --outfile <outfile> 

import sys
import csv
import json
import argparse
import os
from collections import defaultdict
from tqdm import tqdm
import pathogenprofiler as pp
import tbprofiler as tbprofiler

def invert_dict(d):
    inverse = dict()
    for key in d:
        # Go through the list that is saved in the dict:
        for item in d[key]:
            # Check if in the inverted dict the key exists
            if item not in inverse:
                # If not create a new list
                inverse[item] = [key]
            else:
                inverse[item].append(key)
    return inverse

def main(args):

    # tbprofiler_results_location = 'tbprofiler_pakistan_results/'
    metadata = args.metadata
    tbprofiler_results_location = args.tbp_results
    clusters_file = args.clusters_file
    outfile = args.outfile

    # Read in metadata
    # metadata = 'metadata/pakistan_metadata.csv'
    meta_reader = csv.DictReader(open(metadata))
    meta_dict = {}
    for row in meta_reader:
        # Make the id the key, but also recapitulate the id in the key-values by including everything
        meta_dict[row['id_year']] = row

    # Read in clusters file
    # Not sure why this one is read in differently to the metadata
    # {cluster number : [samples list]}
    # clusters_file = 'metadata/PAKISTAN_ALL.clusters.csv'
    clusters_reader = csv.reader(open(clusters_file))
    clusters_dict = defaultdict(list)
    for row in clusters_reader:
        clusters_dict[row[0]].append(row[1])

    # Get first dict in meta_dict
    first_dict = next(iter(meta_dict.values()))

    # Get drug-resistance test names (subset of keys in first_val)
    drug_tests = [x for x in first_dict.keys() if '_test' in x]
    tbprofiler_keys = [x for x in first_dict.keys() if '_tbp' in x]

    # Define the false positive/false negative labels interested in
    outcomes = ["pheno_sens; geno_res", "pheno_res; geno_sens"]
    # Define new dict
    filter_dict = defaultdict(list)
    # Pull together the keys to subset the metadata dictonary
    keys_needed = ['id_year', 'wgs_id'] + drug_tests + tbprofiler_keys
    # Loop through the clusters
    for clust in clusters_dict:
        # Loop through the samples belonging to each cluster
        for samp in clusters_dict[clust]:
            # Loop through the DR tests
            for test in drug_tests:
                # If a sample has a false positive/false negative for any drug, append the subset of the metadata dict to the new dict
                if meta_dict[samp][test] in outcomes:
                    filter_dict[clust].append({key:meta_dict[samp][key] for key in keys_needed})
                    # Break the loop because only need once per sample
                    break

    # Once the samples with false positive/false negative (per cluster) have been found, find the variants that are non-drug resistance associated
    other_variants_dict = {}
    # Loop through the dictionary just created (by cluster)
    for clust in filter_dict:
        # Create empty list per cluster
        other_variants_dict[clust] = []
        # Loop through the dictionaries (one per sample) belonging to each cluster
        for samp_dict in filter_dict[clust]:
            # Get the id and open the json file containing the 'other variants' for it
            id = samp_dict['wgs_id']
            json_file = tbprofiler_results_location + id + '.results.json'
            tbp_result = json.load(open(json_file))
            # Loop over the other_variants dictionaries
            for variant in tbp_result['other_variants']:
                # Exclude synonymous
                if variant['type'] != 'synonymous':
                    # Put it all together
                    other_variants_dict[clust].append({'cluster': clust, 'wgs_id': variant['sample'], 'gene': variant['gene'], 'genome_pos': variant['genome_pos'], 'type': variant['type'],
                    'change': variant['change'], 'nucleotide_change': variant['nucleotide_change'],'locus_tag': variant['locus_tag']})

    # Save a tab-sep text file

    # Define headers from the first dict
    fieldnames = tuple(next(iter(other_variants_dict.values()))[0].keys())

    with open(outfile, 'w') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        # Loop over the dictionaries, appending each dictionary as a row in the file
        for clust in other_variants_dict:
            writer.writerows(other_variants_dict[clust])

parser = argparse.ArgumentParser(description='tbprofiler script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
# parser.add_argument('--samples',type=str,help='File with samples')
parser.add_argument('--metadata',type=str,help='metadata file')
parser.add_argument('--clusters-file',type=str,help='clusters file from ')
parser.add_argument('--tbp-results', default="results/",type=str,help='tbprofiler results directory (json files)')
# parser.add_argument('--db',default="tbdb",type=str,help='Database name')
parser.add_argument('--outfile',default="other_variants.txt",type=str,help='name of output file')
# parser.add_argument('--suffix',default=".results.json",type=str,help='File suffix')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
