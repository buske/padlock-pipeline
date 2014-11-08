#!/usr/bin/env bash

set -eu
set -o pipefail

# Kmer size for repeat identification
kmer=40
# Maximum kmer mappings across the genome
kmer_max_n=5
# Maximum bp overlap between a probe and repetitive regions
max_overlap=50
window=$(python -c "print(${kmer}/2)")

region=region.bed
candidates=probes.bed
mappability="wgEncodeCrgMapabilityAlign${kmer}mer.region.wig"
repeats="${kmer}mer.max${kmer_max_n}.bed"
overlap="probes.max${kmer_max_n}.overlap"
probes="probes.${kmer}mer.max${kmer_max_n}.${max_overlap}bp.txt"
bedgraph=$(basename "${probes}" .txt).bedGraph

# Select region of interest
echo "chr2	136540000	136650000" > "${region}"

# Generate probe file
if [[ ! -e "${candidates}" ]]; then
    cat <<EOF
Error: missing file ${candidates}

Generate candidate probes FASTA for target region using fa2oligo.py,
align to genome using BWA to remove non-unique mapping reads,
and convert to BED.
EOF
    exit 1
fi

# Get mappability track
if [[ ! -e "${mappability}" ]]; then
    cat <<EOF
Error: missing file ${mappability}
Download wgEncodeCrgMapabilityAlign${kmer}mer track from UCSC Table Browser
as 'data points', restricted to region of interest.
EOF
    exit 1
fi

# Identify repetitive regions, for given thresholds
egrep -v "^(track|#)" "${mappability}" \
    | awk -v max_n="${kmer_max_n}" '$4 < (1/max_n)' \
    | cut -f -3 \
    | awk -v n="${window}" 'BEGIN{OFS="\t"} {$2-=n; $3+=n; print $0}' \
    | bedtools merge \
    > "${repeats}"

# Measure coverage with calculated probes
bedtools coverage -b "${candidates}" -a "${repeats}" > "${overlap}"

# Identify non-repetitive probes
awk -v max_bp="${max_overlap}" '$6 < max_bp' "${overlap}" \
    | sort -k 2,2n \
    | cut -f -4 \
    > "${probes}"

# Make bedGraph file for Genome Browser upload
name="${kmer}mer max${kmer_max_n} ${max_overlap}bp"
echo "track type=bedGraph name='${name}' description='${name}'" > "${bedgraph}"
cut -f -3 "${probes}" \
    | bedtools coverage -d -b "${region}" -a stdin \
    | awk 'BEGIN{OFS="\t"}{print $1, $2+$4-1, $2+$4, $5}' \
    >> "${bedgraph}"
