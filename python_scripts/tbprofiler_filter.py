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

    # -----
    # ARGS
    # -----

    # tbprofiler_results_location = 'tbprofiler_pakistan_results/'
    metadata_file = args.metadata_file
    id_key = args.id_key
    tbprofiler_results_location = args.tbp_results
    outfile = args.outfile
    db = args.db

    # -------------
    # READ IN DATA
    # -------------

    # Read in metadata
    meta_reader = csv.DictReader(open(metadata_file))
    meta_dict = {}
    for row in meta_reader:
        # Make the id the key, but also recapitulate the id in the key-values by including everything
        meta_dict[row[id_key]] = row

    # Read in locus-drug resistance associations
    bed_file = "%s/share/tbprofiler/%s.bed" % (sys.base_prefix, db)
    locus_tag2drugs = tbprofiler.get_lt2drugs(bed_file)

    # Get list of files in tbprofiler results directory
    tbprofiler_results_files = os.listdir(tbprofiler_results_location)

    # --------
    # WRANGLE
    # --------

    samples = list(meta_dict.keys())

    # ----------------
    # DR VARIANTS
    # ----------------

    dr_variants_dict = {}
    for json_file in tbprofiler_results_files:
        id = ''.join(json_file.split(".")[:-2])
        if id in samples:
            # Create empty list per id
            dr_variants_dict[id] = []
            json_file = tbprofiler_results_location + json_file
            tbp_result = json.load(open(json_file))
            # Loop over the other_variants dictionaries
            for variant in tbp_result['dr_variants']:
                # print("VARIANT: ", variant['sample'])
                # Exclude synonymous
                if variant['type'] != 'synonymous':
                    # Put it all together. Left join locus/gene drug resistance associations from locus_tag2drugs table
                    empty_str = ""
                    variant.setdefault("gene", empty_str)
                    variant.setdefault("genome_pos", empty_str)
                    variant.setdefault("type", empty_str)
                    variant.setdefault("change", empty_str)
                    variant.setdefault("nucleotide_change", empty_str)
                    variant.setdefault("locus_tag", empty_str)
                    locus_tag2drugs.setdefault(variant['locus_tag'], empty_str)
                    dr_variants_dict[id].append({'wgs_id': id, 'gene': variant['gene'], 'genome_pos': variant['genome_pos'], 'type': variant['type'],
                    'change': variant['change'], 'nucleotide_change': variant['nucleotide_change'], 'locus_tag': variant['locus_tag'], 'locus_tag_drugs': locus_tag2drugs[variant['locus_tag']]})
    
    # Save a tab-sep text file

    # Define headers from the first dict
    fieldnames = tuple(next(iter(dr_variants_dict.values()))[0].keys())

    with open(outfile + '.dr.txt', 'w') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        # Loop over the dictionaries, appending each dictionary as a row in the file
        for id in dr_variants_dict:
            writer.writerows(dr_variants_dict[id])


    # ----------------
    # OTHER VARIANTS
    # ----------------

    other_variants_dict = {}
    for json_file in tbprofiler_results_files:
        id = ''.join(json_file.split(".")[:-2])
        if id in samples:
            # Create empty list per id
            other_variants_dict[id] = []
            json_file = tbprofiler_results_location + json_file
            tbp_result = json.load(open(json_file))
            # Loop over the other_variants dictionaries
            for variant in tbp_result['other_variants']:
                # print("VARIANT: ", variant['sample'])
                # Exclude synonymous
                if variant['type'] != 'synonymous':
                    # Put it all together. Left join locus/gene drug resistance associations from locus_tag2drugs table
                    empty_str = ""
                    variant.setdefault("gene", empty_str)
                    variant.setdefault("genome_pos", empty_str)
                    variant.setdefault("type", empty_str)
                    variant.setdefault("change", empty_str)
                    variant.setdefault("nucleotide_change", empty_str)
                    variant.setdefault("locus_tag", empty_str)
                    locus_tag2drugs.setdefault(variant['locus_tag'], empty_str)
                    # other_variants_dict[clust].append({'wgs_id': id, 'gene': variant['gene'], 'genome_pos': variant['genome_pos'], 'type': variant['type'],
                    # 'change': variant['change'], 'nucleotide_change': variant['nucleotide_change'], 'locus_tag': variant['locus_tag'], 'locus_tag_drugs': locus_tag2drugs[variant['locus_tag']]})
                    other_variants_dict[id].append({'wgs_id': id, 'gene': variant['gene'], 'genome_pos': variant['genome_pos'], 'type': variant['type'],
                    'change': variant['change'], 'nucleotide_change': variant['nucleotide_change'], 'locus_tag': variant['locus_tag'], 'locus_tag_drugs': locus_tag2drugs[variant['locus_tag']]})
    
    # Save a tab-sep text file

    # Define headers from the first dict
    fieldnames = tuple(next(iter(other_variants_dict.values()))[0].keys())

    with open(outfile + '.other.txt', 'w') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        # Loop over the dictionaries, appending each dictionary as a row in the file
        for id in other_variants_dict:
            writer.writerows(other_variants_dict[id])

parser = argparse.ArgumentParser(description='tbprofiler script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--metadata-file', default = '', type = str, help = 'metadata file name with column of sample ids - only processes samples in this file')
parser.add_argument('--id-key', default = '', type = str, help = 'column name in metadata file with sample ids')
parser.add_argument('--db',default="tbdb",type=str,help='prefix to bed file for locus-DR associations')
parser.add_argument('--tbp-results', default="results/",type=str,help='tbprofiler results directory (json files)')
parser.add_argument('--outfile',default="variants",type=str,help='name of output file')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
