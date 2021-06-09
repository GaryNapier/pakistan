#! /usr/bin/env python

# tbprofiler_filter.py --db ../tbdb  <metadata>

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

# def main(args):


tbprofiler_results_location = 'tbprofiler_pakistan_results/'

db = 'tbdb'
bed_file = "%s/share/tbprofiler/%s.bed" % (sys.base_prefix, db)
# bed_file = "%s/share/tbprofiler/%s.bed" % (sys.base_prefix,args.db)
locus_tag2drugs = tbprofiler.get_lt2drugs(bed_file)
# Invert dictionary to make searching by drug easier
locus_tag2drugs_inv = invert_dict(locus_tag2drugs)

# Read in metadata
metadata = 'metadata/pakistan_metadata.csv'
# metadata = args.metadata

meta_reader = csv.DictReader(open(metadata))
# meta_dict = defaultdict(list)
meta_dict = {}
for row in meta_reader:
    meta_dict[row['id_year']] = row

# args.clusters_file
clusters_file = 'metadata/PAKISTAN_ALL.clusters.csv'
clusters_reader = csv.reader(open(clusters_file))
clusters_dict = defaultdict(list)
for row in clusters_reader:
    clusters_dict[row[0]].append(row[1])
clusters_dict_inv = invert_dict(clusters_dict)

# Get first dict in meta_dict
first_dict = next(iter(meta_dict.values()))

# Get drug-resistance test names (subset of keys in first_val)
drug_tests = [x for x in first_dict.keys() if '_test' in x]
tbprofiler_keys = [x for x in first_dict.keys() if '_tbp' in x]

outcomes = ["pheno_sens; geno_res", "pheno_res; geno_sens"]
filter_dict = defaultdict(list)
keys_needed = ['id_year', 'wgs_id'] + drug_tests + tbprofiler_keys
for clust in clusters_dict:
    # print("CLUST", clust)
    for samp in clusters_dict[clust]:
        # print("SAMP", samp)
        for test in drug_tests:
            # print("TEST", test)
            if meta_dict[samp][test] in outcomes:
                # print(clust, samp, meta_dict[samp]['wgs_id'], test, meta_dict[samp][test])
                filter_dict[clust].append({key:meta_dict[samp][key] for key in keys_needed})
                break


other_variants_dict = {}
for dict_list in filter_dict:
    # print(filter_dict[dict_list])
    other_variants_dict[dict_list] = []
    for samp_dict in filter_dict[dict_list]:
        # print(samp_dict)
        id = samp_dict['wgs_id']
        print(id)
        json_file = tbprofiler_results_location + id + '.results.json'
        tbp_result = json.load(open(json_file))
        for variant in tbp_result['other_variants']:
            if variant['type'] != 'synonymous':
                other_variants_dict[dict_list].append({'wgs_id': variant['sample'], 'gene': variant['gene'], 'genome_pos': variant['genome_pos'], 'type': variant['type'],
                'change': variant['change'], 'nucleotide_change': variant['nucleotide_change'],'locus_tag': variant['locus_tag']})
other_variants_dict

other_variants_dict['3']




    if args.samples:
        samples = [x.rstrip() for x in open(args.samples).readlines()]
    else:
        samples = [x.replace(args.suffix,"") for x in os.listdir(args.dir) if x[-len(args.suffix):]==args.suffix]

    for s in tqdm(samples):
        data = json.load(open(pp.filecheck("%s/%s%s" % (args.dir,s,args.suffix))))


parser = argparse.ArgumentParser(description='tbprofiler script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
# parser.add_argument('--samples',type=str,help='File with samples')
parser.add_argument('--metadata',type=str,help='Metadata file')
parser.add_argument('--dir',default="results/",type=str,help='Directory containing results')
parser.add_argument('--db',default="tbdb",type=str,help='Database name')
parser.add_argument('--suffix',default=".results.json",type=str,help='File suffix')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
