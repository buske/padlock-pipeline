#!/usr/bin/env python

"""
Combine Bismark methylation coverage files, either by aggregating coverage 
or by appending columns
"""

import os
import sys
import re

from collections import defaultdict


def read_file(filename, region=None, region_re=re.compile('^\s*(\w+)\s*:\s*(\d+)\s*-\s*(\d+)\s*$')):
    if region:
        region_match = region_re.match(region)
        assert region_match, 'Invalid region: {}'.format(region)
        region_chrom, region_start, region_end = region_match.groups()
        region_start = int(region_start)
        region_end = int(region_end)

    counts = {}
    with open(filename) as ifp:
        for line in ifp:
            line = line.rstrip('\n')
            if not line or line.startswith('#'): continue

            tokens = line.split('\t')
            if len(tokens) == 5:
                chrom, start, cpg, n_mc, n_c = tokens
            else:
                raise NotImplementedError()

            start = int(start)
            n_mc = int(n_mc)
            n_c = int(n_c)

            if region and (chrom != region_chrom or
                           start < region_start or
                           start > region_end):
                continue

            locus = (chrom, start, cpg)
            assert locus not in counts
            counts[locus] = (n_mc, n_c)

    return counts

def read_many_files(filenames, aggregate=False, region=None):
    all_loci = set()
    all_counts = []
    for filename in filenames:
        counts = read_file(filename, region=region)
        all_loci.update(counts)
        all_counts.append(counts)
        
    merged_counts = {}
    for locus in all_loci:
        if aggregate:
            row = [0, 0]
        else:
            row = []

        for counts in all_counts:
            count = counts.get(locus, (0, 0))
            if aggregate:
                row[0] += count[0]
                row[1] += count[1]
            else:
                row.extend(count)

        merged_counts[locus] = row

    return merged_counts
    
def script(filenames, aggregate=False, add_header=False, region=None):
    counts = read_many_files(filenames, aggregate=aggregate, region=region)
    loci = sorted(counts)

    if add_header:
        colnames = ['chrom', 'pos', 'cpg']
        if aggregate:
            colnames.extend(['methyl', 'unmethyl'])
        else:
            for filename in filenames:
                base = os.path.dirname(filename)
                colnames.extend([base + '_methyl', base + '_unmethyl'])
        print('#{}'.format('\t'.join(colnames)))
    for locus in loci:
        row = counts[locus]
        line = list(locus)
        line.extend(row)
        print('\t'.join(map(str, line)))

def parse_args(args):
    from argparse import ArgumentParser
    description = __doc__.strip()
    
    parser = ArgumentParser(description=description)
    parser.add_argument('filenames', metavar='BISMARK.COV', nargs='+')
    parser.add_argument('--add-header', default=False, action='store_true')
    parser.add_argument('--region', help='Filter to only sites within given'
                        ' chr:start-end region.')
    parser.add_argument('--aggregate', default=False, action='store_true',
                        help='Add together coverage counts, rather than'
                        ' appending columns')

    return parser.parse_args(args)

def main(args=sys.argv[1:]):
    args = parse_args(args)
    script(**vars(args))

if __name__ == '__main__':
    sys.exit(main())
