#!/usr/bin/env python

"""
Given a single-sequence fasta file, report the positions of all CpGs.
"""

__author__ = 'Orion Buske (buske@cs.toronto.edu)'

import os
import sys
import logging
import re

def maybe_stdin_open(filename):
    if filename == '-':
        return sys.stdin
    else:
        return open(filename)

def read_fasta(filename):
    parts = []
    with maybe_stdin_open(filename) as ifp:
        header = ifp.readline()
        assert header.startswith('>')
        seq_name = header.lstrip('>').strip()
        logging.info('Loading sequence : {}'.format(seq_name))
        for line in ifp:
            line = line.rstrip('\n')
            assert not line.startswith('>'), \
                'Found more than one sequence in FASTA'
            parts.append(line)

    return seq_name, ''.join(parts)

def script(fasta_filename, min_position, max_position):
    seq_name, seq = read_fasta(fasta_filename)
    logging.info('Found sequence {!r}: {}bp'.format(seq_name, len(seq)))

    if min_position is None:
        min_position = 1

    if max_position is None:
        max_position = len(seq)

    logging.info('Finding CpGs within: [{}, {}]'.format(min_position, max_position))


    cpg_re = re.compile(r'CG', re.IGNORECASE)
    for match in cpg_re.finditer(seq, min_position - 1, max_position):
        # match.start is index of C, convert to 1-based position
        cpg_position = match.start() + 1
        print('{}\t{}'.format(seq_name, cpg_position))

def parse_args(args):
    from argparse import ArgumentParser
    description = __doc__.strip()
    
    parser = ArgumentParser(description=description)
    parser.add_argument('fasta_filename', metavar='FASTA',
                        help='Single-sequence fasta ("-" to read from stdin)')
    parser.add_argument('--min', dest='min_position', type=int,
                        help='Minimum position to report')
    parser.add_argument('--max', dest='max_position', type=int,
                        help='Maximum position to report')
    parser.add_argument('--log', default='INFO',
                        help='Logging level (e.g. DEBUG, INFO, WARNING)')

    return parser.parse_args(args)

def main(args=sys.argv[1:]):
    args = vars(parse_args(args))
    logging.basicConfig(level=args.pop('log'))
    script(**args)

if __name__ == '__main__':
    sys.exit(main())
