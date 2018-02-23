#!/bin/python
"""
selective parse for UniProt::IRGSP mappings
"""

import argparse
import re

# process args
parser = argparse.ArgumentParser(description='Selective parse for UniProt::IRGSP mappings')

# input settings
parser.add_argument('-i', '--input_path', help='2-col tab file with uniprot and rice genes')
args = parser.parse_args()

# regex for gene name matching
valid_IRGSP = re.compile("OS..G.......")
valid_MSU = re.compile("LOC_OS..G.....")

count = 0

if args.input_path:
    # iter input file
    INF = open(args.input_path)
    for line in INF:
        cols = line.rstrip().split('\t')
        uniprot = cols[0]
        # grab col 2 and split on space
        os_gene_names = cols[1].lstrip(' ').upper().split(' ')
        found = 0
        # iter list and break on first IRGSP match: os..g.......
        for os_gene_name in os_gene_names:
            if valid_IRGSP.match(os_gene_name):
                found = 1
                count += 1
                print(uniprot + "\t" + os_gene_name)
                break
        # if no IRGSP match, re-iter and break on first MSU match
        if not found:
            for os_gene_name in os_gene_names:
                if valid_MSU.match(os_gene_name):
                    found = 1
                    count += 1
                    print(uniprot + "\t" + os_gene_name)
                    break
        if not found:
            print("NO MAP")
    INF.close()

#print("map #: " + str(count))
