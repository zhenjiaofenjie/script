#!/usr/bin/env python
# -*- coding: utf-8 -*-

from Bio import SeqIO
import argparse

#read the arguments
parser = argparse.ArgumentParser(description='Read GenBank file and get the new and old locus tag.')
parser.add_argument('input', type=str, help='Input file name')
parser.add_argument('output', type=str, help='Output file name')
args = parser.parse_args()

with open(args.output, 'w') as out:
    for record in SeqIO.parse(args.input, "genbank"):
        for ft in record.features:
            if ft.type == 'CDS':
                if 'old_locus_tag' in ft.qualifiers: #Some genes don't exist in older version
                    out.write('{}\t{}\n'.format(ft.qualifiers['locus_tag'][0], ft.qualifiers['old_locus_tag'][0]))
                else:
                    out.write('{}\tNA\n'.format(ft.qualifiers['locus_tag'][0]))
