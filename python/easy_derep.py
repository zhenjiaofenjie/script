#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse, textwrap
import csv

class Table(object):
    def __init__(self):
        self.table = []

    def ReadTable(self, inputtable, column):
        keylist = []
        with open(inputtable, 'r') as i:
            reader = csv.reader(i, dialect='excel-tab')
            for linelist in reader:
                key = []
                for col in column:
                    key.append(linelist[col-1]) #make the index user-friendly
                if key in keylist:
                    continue
                self.table.append(linelist)
                keylist.append(key)

    def PrintTable(self, output):
        with open(output, 'w') as o:
            wr = csv.writer(o, delimiter= '\t')
            wr.writerows(self.table)

if __name__ == '__main__':
    #read the arguments
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent('''
        Delete duplicate rows based on selected columns
        '''))
    parser.add_argument('-i', '--input', type=str, help='Input table', action='store', required=True)
    parser.add_argument('-c', '--columns', type=int, help='The index of the key column in the table', action='append', required=True)
    parser.add_argument('-o', '--output', type=str, help='Output file name', required=True)
    args = parser.parse_args()

    table = Table()
    table.ReadTable(args.input, args.columns)
    table.PrintTable(args.output)
