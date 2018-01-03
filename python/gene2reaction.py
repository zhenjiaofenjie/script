#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse, textwrap
import csv, re

class Reference(object):
    def __init__(self):
        self.refdist = {}

    def LoadReference(self, rfile, rc, gene):
        with open(rfile, 'r') as r:
            reader = csv.reader(r, dialect='excel-tab')
            for row in reader:
                reaction = row[rc-1]
                genelist = row[gene-1]
                genelist = re.sub('[,\(\)(or)(and)]*', '', genelist)
                genelist = re.split('\s+', genelist)
                for i in genelist:
                    if i in self.refdist:
                        self.refdist[i].append(reaction)
                    else:
                        self.refdist[i]=[reaction]

    def WriteReference(self, output):
        with open(output, 'w') as o:
            for k, v in self.refdist.iteritems():
                o.write('{}\t{}\n'.format(k, ','.join(v)))

if __name__ == '__main__':
    #read the arguments
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent('''
        Read csv file listing the genes for each reaction.
        Output the corresponding reactions for each gene.
        '''))
    parser.add_argument('-i', '--input', type=str, help='Reaction to gene file. Can contain [,|and|or] in gene list.', action='store', required=True)
    parser.add_argument('-r', '--rc', type=int, help='The index of the reaction column in the table', required=True)
    parser.add_argument('-g', '--gene', type=int, help='The index of the gene list column in the table', required=True)
    parser.add_argument('-o', '--output', type=str, help='Output file name', required=True)
    args = parser.parse_args()

    refdist = Reference()
    refdist.LoadReference(args.input, args.rc, args.gene)
    refdist.WriteReference(args.output)
