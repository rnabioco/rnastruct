#!/usr/bin/env bash
#BSUB -J RNAstruct 
#BSUB -o log/snakemake_%J.out
#BSUB -e log/snakemake_%J.err
#BSUB -R "select[mem>4] rusage[mem=4]"
#BSUB -q rna
set -o nounset -o pipefail -o errexit -x


args=' -q rna -o {log}.out -e {log}.err -J {params.job_name} -R "{params.memory}" -n {threads}'

. /usr/share/Modules/init/bash
module load modules modules-init modules-python
module load anaconda
unset PYTHONPATH
source activate py37

snakemake \
    --drmaa "$args" \
    --snakefile Snakefile \
    --jobs 160 \
    --resources all_threads=160 \
    --latency-wait 300 \
    --rerun-incomplete \
    --configfile config.yaml
    #--until deeptools_norm_factors 
