
fa="$HOME/Projects/shared_dbases/genomes/human/chr_appended_grch37/Homo_sapiens.GRCh37.dna.primary_assembly_chrappended.fa"
tbl_wt="$HOME/Projects/rnastruct/data/mismatches/june2019/filtered_tables/50/WT/pileup_table.sorted.tsv.bgz"
tbl_c4="$HOME/Projects/rnastruct/data/mismatches/june2019/filtered_tables/50/C4/pileup_table.sorted.tsv.bgz"
bed="examples.bed"

extract_transcript.py \
    -t $tbl_wt \
    -g $bed \
    -f $fa \
    -m 24 \
    -e 25 \
    -o wt_map

extract_transcript.py \
    -t $tbl_c4 \
    -g $bed \
    -f $fa \
    -m 24 \
    -e 25 \
    -o c4_map
