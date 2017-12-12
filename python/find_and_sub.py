#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse, textwrap
import csv

def LoadReference(ref):
    refdist = {}
    for rfile in ref:
        with open(rfile, 'r') as r:
            reader = csv.reader(r, dialect='excel-tab')
            for row in reader:
                refdist[row[0]] = row[1]
    return refdist

def FindAndSub(input, refdict, output):
    with open(input, 'r') as f:
        with open(output, 'w') as o:
            reader = csv.reader(f, dialect='excel-tab')
            for row in reader:
                newrow = []
                for cell in row:
                    newcell = []
                    for i in cell.split(','):
                        if i in refdict:
                            i = refdict[i]
                        elif i is not '-':
                            print('Warning: {} is not in reference.'.format(i))
                        newcell.append(i)
                    newrow.append(','.join(newcell))
                o.write('\t'.join(newrow))
                o.write('\n')

if __name__ == '__main__':
    #read the arguments
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent('''
        Find and substitution.
        Multiple reference files is allowed.
        '''))
    parser.add_argument('-i', type=str, help='Input file name', required=True, metavar='INPUT')
    parser.add_argument('-r', type=str, help='Reference file name[s]', action='append', required=True, metavar='RF')
    parser.add_argument('-o', '--output', type=str, help='Output file name', required=True)
    args = parser.parse_args()

    refdist = LoadReference(args.r)
    FindAndSub(args.i, refdist, args.output)
