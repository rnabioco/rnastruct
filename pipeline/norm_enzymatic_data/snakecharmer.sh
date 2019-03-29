#!/usr/bin/env bash
#BSUB -J RNAstruct 
#BSUB -o log/snakemake_%J.out
#BSUB -e log/snakemake_%J.err
#BSUB -R "select[mem>4] rusage[mem=4]"
#BSUB -q rna
set -o nounset -o pipefail -o errexit -x


args=' -q rna -o {log}.out -e {log}.err -J {params.job_name} -R "{params.memory}" -n {threads}'

snakemake \
    --drmaa "$args" \
    --snakefile Snakefile \
    --jobs 36 \
    --resources all_threads=36 \
    --latency-wait 300 \
    --rerun-incomplete \
    --configfile config.yaml 
