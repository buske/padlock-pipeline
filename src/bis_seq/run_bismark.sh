#!/usr/bin/env bash

set -eu
set -o pipefail


##### Global parameters #####
suffix=.fastq.gz
thisdir="$(cd "$(dirname $0)"; pwd)"
srcdir="$(pwd)/src"
refdir="$(pwd)/sequences"
trimmomatic_jar="$(pwd)/lib/Trimmomatic-0.32/trimmomatic-0.32.jar"
trimmomatic="java -Xmx2048m -jar ${trimmomatic_jar}"
bismarkdir="$(pwd)/lib/bismark_v0.12.3"
bowtiedir="$(pwd)/lib/bowtie-1.1.1"
adapters="${refdir}/adapters.fa"

# DISPATCH: method for running SGE-compatible subscripts
DISPATCH="${DISPATCH:-qsub}"

##### Parse command-line arguments and ensure dependencies exist #####
function usage {
    cat <<EOF
Usage: $0 ASSEMBLY DIR

Script must be run from the root of the padlock-pipeline repo.

ASSEMBLY: the assembly to which the sample will be aligned. This must be a
  subdirectory of ${refdir}/bismark, already prepared with 
  "bismark_genome_preparation".

DIR: a sample directory, containing read files with names of the format 
  "<SAMPLE>_L001_R1_001${suffix}". 

- ${refdir}/bismark/phiX must already exist and be prepared with bismark.
- All reads from the same strand will be merged across lanes and file parts.
- The script creates SGE-compatible bash scripts within the sample directory.
- Depending on the value of the DISPATCH environment variable, scripts will
  be run over SGE, run on the local machine, or not run:
    'qsub': scripts will be automatically dispatched 
    'bash': scripts will be in serial on the local machine
    'echo': scripts will not be run (their paths will be echoed
EOF
    exit 1
}

function require_files {
    while [[ $# -gt 0 ]]; do
	file="$1"
	shift
	if [[ ! -e "${file}" ]]; then
	    echo "Unable to find required file: ${file}" >&2
	    exit 1
	fi
    done
}

function require_dirs {
    while [[ $# -gt 0 ]]; do
	dir="$1"
	shift
	if [[ ! -d "${dir}" ]]; then
	    echo "Unable to find required directory: ${dir}" >&2
	    exit 1
	fi
    done
}

require_files "${trimmomatic_jar}" "${adapters}"
require_dirs "${srcdir}" "${refdir}" "${bismarkdir}" "${bowtiedir}"

if [[ DISPATCH = "qsub" ]]; then
    which qsub > /dev/null
elif [[ "${DISPATCH}" = "echo" ]]; then
    require_dir "${TMPDIR}"
fi

if [[ $# -eq 0 ]]; then
    usage
elif [[ $# -ne 2 ]]; then
    echo -e "Error: unexpected number of arguments\n" >&2
    usage
fi


export PATH="${bismarkdir}:${bowtiedir}:${PATH}"


##### Define functions #####

function trim_sample {
    local lanes=( $(cd "${dir}" && ls *"${suffix}" | perl -pe 's/.*(L\d+).*/\1/' | sort -u) )

    if [[ "${#lanes[@]}" -lt 1 ]]; then
	echo "Unable to find lanes for ${name}" >&2
	return 1
    fi

    if [[ ! -s "${dir}/${name}.R1P.fq.gz" ]]; then
	local script="${dir}/pipeline-trim.sh"
	cat > "${script}" <<EOF
#!/usr/bin/env bash

#$ -pe parallel 2
#$ -l h_vmem=6G
#$ -V
#$ -S /bin/bash
#$ -N "${name}.trim"
#$ -e "${logdir}"
#$ -o "${logdir}"
#$ -l hostname="supa*"

set -eu

cd "\${TMPDIR}"


# Paired-end trimming
for lane in ${lanes[@]}; do
  r1="\$(ls "${dir}"/*_\${lane}_R1_*${suffix})"
  r2="\$(ls "${dir}"/*_\${lane}_R2_*${suffix})"
  test -s "\${r1}"
  test -s "\${r2}"

  ${trimmomatic} PE -threads 2 -phred33 "\${r1}" "\${r2}" \
    "${name}.\${lane}.R1P.fq.gz" "${name}.\${lane}.R1U.fq.gz" "${name}.\${lane}.R2P.fq.gz" "${name}.\${lane}.R2U.fq.gz" \
    ILLUMINACLIP:"${adapters}":2:25:7:1 LEADING:10 TRAILING:10 SLIDINGWINDOW:3:15 MINLEN:75
  cat "${name}.\${lane}.R1P.fq.gz" >> "${name}.R1P.fq.gz"
  cat "${name}.\${lane}.R2P.fq.gz" >> "${name}.R2P.fq.gz"
  cat "${name}.\${lane}.R1U.fq.gz" >> "${name}.R1U.fq.gz"
  cat "${name}.\${lane}.R2U.fq.gz" >> "${name}.R2U.fq.gz"
done

rm -f "${name}".*.*.fq.gz

mv *.fq.gz "${dir}/"
EOF

	${DISPATCH} "${script}"
    fi
}

function unspike {
    local assembly=phiX
    local ref="${refdir}/bismark/$assembly"

    if [[ ! -e "${dir}/${name}.${assembly}.paired_pe.sam" || ! -e "${dir}/${name}.${assembly}.unpaired_r1.sam" || ! -e "${dir}/${name}.${assembly}.unpaired_r2.sam" ]]; then
	local script="${dir}/pipeline-align-${assembly}.sh"
	cat > "${script}" <<EOF
#!/usr/bin/env bash

#$ -l h_vmem=16G
#$ -V
#$ -S /bin/bash
#$ -hold_jid "${name}.trim"
#$ -N "${name}.unspike"
#$ -e "${logdir}"
#$ -o "${logdir}"
#$ -l hostname="supa*"

set -eu

cd "\${TMPDIR}"

if [[ ! -e "${dir}/${name}.${assembly}.paired_pe.sam" ]]; then
  "${bismarkdir}/bismark" --path_to_bowtie "${bowtiedir}/" --output_dir . --temp_dir . --basename "${name}.${assembly}.paired" --unmapped -n 1 --non_directional -X 1000 "${ref}" -1 "${dir}/${name}.R1P.fq.gz" -2 "${dir}/${name}.R2P.fq.gz"
fi

for r in 1 2; do
  if [[ ! -e "${dir}/${name}.${assembly}.unpaired_r\${r}.sam" ]]; then
    "${bismarkdir}/bismark" --path_to_bowtie "${bowtiedir}/" --output_dir . --temp_dir . --basename "${name}.${assembly}.unpaired_r\${r}" --unmapped -n 1 --non_directional "${ref}" "${dir}/${name}.R\${r}U.fq.gz"
  fi
done

mv "${name}.${assembly}".* "${dir}/"
EOF
	${DISPATCH} "${script}"
    fi
}


function align {
    if [[ ! -e "${dir}/${name}.${assembly}.paired_pe.sam" || ! -e "${dir}/${name}.${assembly}.unpaired_r1.sam" || ! -e "${dir}/${name}.${assembly}.unpaired_r2.sam" ]]; then
	local script="${dir}/pipeline-align-${assembly}.sh"
	cat > "${script}" <<EOF
#!/usr/bin/env bash

#$ -l h_vmem=16G
#$ -V
#$ -S /bin/bash
#$ -hold_jid "${name}.unspike"
#$ -N "${name}.align.${assembly}"
#$ -e "${logdir}"
#$ -o "${logdir}"
#$ -l hostname="supa*"

set -eu

cd "\${TMPDIR}"

if [[ -e "${dir}/${name}.phiX.paired_unmapped_reads_1.fq" && -e "${dir}/${name}.phiX.paired_unmapped_reads_2.fq" && ! -e "${dir}/${name}.${assembly}.paired_pe.sam" ]]; then
  "${bismarkdir}/bismark" --path_to_bowtie "${bowtiedir}/" --output_dir . --temp_dir . --basename "${name}.${assembly}.paired" -n 1 --non_directional -X 1000 "${ref}" -1 "${dir}/${name}.phiX.paired_unmapped_reads_1.fq" -2 "${dir}/${name}.phiX.paired_unmapped_reads_2.fq"
fi

for r in r1 r2; do
  if [[ -e "${dir}/${name}.phiX.unpaired_\${r}_unmapped_reads.fq" && ! -e "${dir}/${name}.${assembly}.unpaired_\${r}.sam" ]]; then
    "${bismarkdir}/bismark" --path_to_bowtie "${bowtiedir}/" --output_dir . --temp_dir . --basename "${name}.${assembly}.unpaired_\${r}" -n 1 --non_directional "${ref}" "${dir}/${name}.phiX.unpaired_\${r}_unmapped_reads.fq"
  fi
done

mv "${name}.${assembly}".* "${dir}/"
EOF
	${DISPATCH} "${script}"
    fi
}

function coverage {
    if [[ ! -s "${dir}/${name}.${assembly}.merged.cov" ]]; then
	local script="${dir}/pipeline-coverage-${assembly}.sh"
	cat > "${script}" <<EOF
#!/usr/bin/env bash

#$ -l h_vmem=8G
#$ -V
#$ -S /bin/bash
#$ -hold_jid "${name}.align.${assembly}"
#$ -N "${name}.cov.${assembly}"
#$ -e "${logdir}"
#$ -o "${logdir}"
#$ -l hostname="supa*"

set -eu
set -o pipefail

cd "\${TMPDIR}"

# Generate any missing coverage files
for pairing in paired_pe unpaired_r1 unpaired_r2; do
  if [[ ! -s "${dir}/${name}.${assembly}.\${pairing}.bismark.cov" ]]; then
    if [[ "\${pairing}" = "paired_pe" ]]; then
      extra_args="-p --no_overlap"
    else
      extra_args=""
    fi

    "${bismarkdir}/bismark_methylation_extractor" \${extra_args} --comprehensive --merge_non_CpG --report "${dir}/${name}.${assembly}.\${pairing}.sam"
    "${srcdir}/mybismark2bedGraph" --CX_context --buffer_size 4G -o "${name}.${assembly}.\${pairing}.bedGraph" {,Non_}"CpG_context_${name}.${assembly}.\${pairing}.txt"
    mv -v "${name}.${assembly}.\${pairing}.bedGraph" "${dir}/"
    mv -v "${name}.${assembly}.\${pairing}.bismark.cov" "${dir}/"
  fi
done

# Sum coverages for paired and unpaired reads
"${srcdir}/merge_bismark_coverage.py" --aggregate "${dir}/${name}.${assembly}".{paired_pe,unpaired_r1,unpaired_r2}.bismark.cov > "${name}.${assembly}.merged.cov"

mv -v "${name}.${assembly}.merged.cov" "${dir}/"

EOF
	${DISPATCH} "${script}"
    fi
}

function process_dir {
    # Set some globals that most functions will need to use
    assembly="$1"
    dir="$2"

    name="$(basename "${dir}")"
    name="${name#Sample_}"

    logdir="${dir}/log"
    mkdir -pv "${logdir}"
    
    ref="${refdir}/bismark/${assembly}"

    # Step 1a: trim
    trim_sample

    # Step 1b: remove phiX
    unspike

    # Step 2: align with bismark
    align

    # Step 3: measure methylation coverage
    coverage
}


##### Run program #####
process_dir "$1" "$2"