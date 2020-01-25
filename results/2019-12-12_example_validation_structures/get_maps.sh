
src="$HOME/Projects/rnastruct/src"
fa="$HOME/Projects/shared_dbases/genomes/human/chr_appended_grch37/Homo_sapiens.GRCh37.dna.primary_assembly_chrappended.fa"

tbl_dir="$HOME/Projects/rnastruct/data/mismatches/oct2019/all_libs/filtered_tables_no_cutoff/15"
tbl_wt=$tbl_dir"/WT/pileup_table.sorted.tsv.bgz"
tbl_c4=$tbl_dir"/C4_subset/pileup_table.sorted.tsv.bgz"
bed="examples.bed"

$src/extract_transcript.py \
    -t $tbl_wt \
    -g $bed \
    -f $fa \
    -m 19 \
    -e 20 \
    -o wt_map

$src/extract_transcript.py \
    -t $tbl_c4 \
    -g $bed \
    -f $fa \
    -m 19 \
    -e 20 \
    -o c4_map
