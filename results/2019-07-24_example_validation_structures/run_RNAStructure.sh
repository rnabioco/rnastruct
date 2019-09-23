

in=$1

python ~/Projects/rnastruct/bin/mapToFastaShape.py \
    --uConvert \
    $in

Fold \
    ${in/.map/.fa} \
    ${in/.map/.db} \
    -k \
    -sh ${in/.map/.shape} \
    -sm 1.8 \
    -si -0.6 \
    -md 600 \
    -m 100 
