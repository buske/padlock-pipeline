#!/usr/bin/env python

"""
Given a single-sequence fasta, generate candidate oligonuclotide sequences
around restriction sites in a fasta output.
"""

__author__ = 'Orion Buske (buske@cs.toronto.edu)'

import os
import sys
import logging

def read_fasta(filename):
    parts = []
    with open(filename) as ifp:
        header = ifp.readline()
        assert header.startswith('>')
        seq_name = header.lstrip('>').strip()
        logging.info('Designing probes for sequence: {}'.format(seq_name))
        for line in ifp:
            line = line.rstrip('\n')
            assert not line.startswith('>')
            parts.append(line)

    return ''.join(parts)

def generate_probes(fragment, probe_length, max_probe_distance):
    probes = set()
    if len(fragment) < probe_length:
        return probes

    # Number of steps to tile from each end
    n_tiles = min(len(fragment), max_probe_distance) - probe_length
    logging.debug('Tiling {} probes from each end'.format(n_tiles))

    for offset in range(n_tiles):
        oligo = fragment[offset:offset + probe_length]
        probes.add(oligo)

    # If smaller than max_probe_distance, probes already span full fragment
    if len(fragment) > max_probe_distance:
        # Otherwise, need to capture some probes from other end
        for offset in range(n_tiles):
            start = len(fragment) - offset - probe_length
            end = len(fragment) - offset
            oligo = fragment[start:end]
            probes.add(oligo)

    return probes

def print_probes(probes):
    for i, probe in enumerate(probes):
        print('>probe_{}\n{}'.format(i, probe))

def script(fasta_filename, restriction_site, probe_length, max_probe_distance):
    seq = read_fasta(fasta_filename)
    logging.info('Sequence is {}bp'.format(len(seq)))

    fragments = seq.split(restriction_site)
    logging.info('{} fragments expected after restriction with {}'.format(len(fragments), restriction_site))

    probes = set()
    for fragment in fragments:
        logging.debug('Fragment of length: {}'.format(len(fragment)))
        frag_probes = generate_probes(fragment, probe_length, max_probe_distance)
        logging.debug('Designed {} probes'.format(len(frag_probes)))
        probes.update(frag_probes)
        
    logging.info('Generated {} unique {}bp probe sequences'.format(len(probes), probe_length))
    print_probes(probes)

def parse_args(args):
    from argparse import ArgumentParser
    description = __doc__.strip()
    
    parser = ArgumentParser(description=description)
    parser.add_argument('fasta_filename', metavar='FASTA',
                        help='Single-sequence fasta to design probes for')
    parser.add_argument('--restriction-site', default='GATC',
                        help='Restriction site sequence (default for DpnII)')
    parser.add_argument('--probe-length', default=120,
                        help='The length of each probe oligo')
    parser.add_argument('--max-probe-distance', default=140,
                        help='Maximum distance from far end of probe to'
                        ' nearest restriction site')
    parser.add_argument('--log', default='INFO',
                        help='Logging level (e.g. DEBUG, INFO, WARNING)')

    return parser.parse_args(args)

def main(args=sys.argv[1:]):
    args = vars(parse_args(args))
    logging.basicConfig(level=args.pop('log'))

    script(**args)

if __name__ == '__main__':
    sys.exit(main())
