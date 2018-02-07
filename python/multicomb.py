#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse, textwrap
import csv
import re

def check_args(args):
    if not len(args.l)==len(args.r)==len(args.lci)==len(args.rci):
        quit('Error: -l -r --lci and --rci should be defined the same times!')

    l_list=[]
    for l, r in zip(args.l, args.r):
        if len(l_list) > 0:
            if l not in l_list:
                quit('Error: {} has not been mentioned before, it can not be the left table yet!'.format(l))
        l_list.append(l)
        l_list.append(r)

class TableCombined(object):
    def __init__(self):
        self.table = []
        self.index = {}
        self.lcidict = {}
        self.rcidict = {}

    def PrintTable(self, output):
        with open(output, 'w') as o:
            wr = csv.writer(o, delimiter= '\t')
            wr.writerows(self.table)

def CombineTable(l, r, lci, rci, combined, ignorecase):
    dictid = {}
    with open(r, 'r') as rtable:
        reader = csv.reader(rtable, dialect='excel-tab')
        for linelist in reader:
            if re.match('^#', linelist[0]):
                continue
            if ignorecase:
                dictid[linelist.pop(rci).lower()] = linelist
            else:
                dictid[linelist.pop(rci)] = linelist
    rcols = len(linelist)


    if l in combined.index: #the left table l has been combined in previous circle
        if l+str(lci) in combined.lcidict: #the lci is the previously omitted key column, use previous lci instead
            lci = combined.lcidict[l+str(lci)]
        elif l in combined.rcidict and lci > combined.rcidict[l]:  #the lci is after the previously omitted key column, thus it should minus the omitted column
            lci += combined.index[l] - 1 #index is to show where the l starts in the combined table
        else:  #the lci is before the previously omitted key column
            lci += combined.index[l]
    else:  #the first time to run CombineTable
        with open(l, 'r') as ltable:
            reader = csv.reader(ltable, dialect='excel-tab')
            for linelist in reader:
                if re.match('^#', linelist[0]):
                    continue
                combined.table.append(linelist)
        combined.index[l] = 0

    combined.lcidict[r+str(rci)] = lci #this column has been omitted, could not be used as key in the future
    combined.rcidict[r] = rci #record the index of omitted key column
    temp = []
    for row in combined.table:
        if row[lci] in dictid:
            temp.append(row+dictid[row[lci]])
        elif row[lci].lower() in dictid and ignorecase:
            temp.append(row+dictid[row[lci].lower()])
        else:
            for i in range(rcols):
                row.append('-')
            temp.append(row)

    lcols = len(combined.table[0])
    combined.table = temp
    combined.index[r] = lcols  #record where the r starts in the combined table
    return combined

if __name__ == '__main__':
    #read the arguments
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent('''
        Combine tables according to given columns.
        It combines a pair of left and right table each time,
        according to the combination of "-l LF -r RF --lci LN --rci RN".
        Use "-l LF -r RF --lci LN --rci RN" muliple times to merge all tables together.
        Make sure each LF except the first one has been mentioned in previous combinations.
        '''))
    parser.add_argument('--ignorecase', action='store_true', help='Ignore the case of key columns')
    parser.add_argument('-l', type=str, help='Left table', action='append', required=True, metavar='LF')
    parser.add_argument('-r', type=str, help='Right table', action='append', required=True, metavar='RF')
    parser.add_argument('--lci', type=int, help='The index of the key column in the left table', action='append', metavar='LN', required=True)
    parser.add_argument('--rci', type=int, help='The index of the key column in the right table', action='append', metavar='RN', required=True)
    parser.add_argument('-o', '--output', type=str, help='Output file name', required=True)
    args = parser.parse_args()

    check_args(args)

    combined = TableCombined()
    for l, r, lci, rci in zip(args.l, args.r, args.lci, args.rci): #args.lci and args.rci start from 1
        combined = CombineTable(l, r, lci-1, rci-1, combined, args.ignorecase) #lci and rci start from 0

    combined.PrintTable(args.output)
