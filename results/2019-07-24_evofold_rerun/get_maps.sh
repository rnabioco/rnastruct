#BSUB -n 1
#BSUB -J filter[1-4]
#BSUB -e err_%J.txt
#BSUB -o out_%J.txt
#BSUB -R "rusage[mem=10] select[mem>10] span[hosts=1]"
#BSUB -q rna


DEPTHS=(
10
15
50
100
)

depth=${DEPTHS[$(($LSB_JOBINDEX - 1))]}
fa="/beevol/home/riemondy/Projects/shared_dbases/genomes/human/chr_appended_grch37/Homo_sapiens.GRCh37.dna.primary_assembly_chrappended.fa"
datadir="/beevol/home/riemondy/Projects/rnastruct/data/mismatches/july2019/filtered_tables/"$depth
cutoff_dir="no_slop"
outdir=$cutoff_dir"/"$depth"_cfiles"

mkdir -p $outdir

extract_transcript.py \
    -t $datadir"/WT/pileup_table.sorted.tsv.bgz" \
    -g $cutoff_dir"/"$depth"_cutoff.txt" \
    -f $fa \
    -o $outdir \
    -T "C" \
    -m 22 \
    -e 23 \
    -n

#mkdir $outdir"/no_dms"
#for file in $outdir/*.map
#do echo $file
#    outfile=$(basename $file)
#    awk '{OFS=FS="\t"} {print $1,-999.0,$3,$4}' $file > $outdir"/no_dms/"$outfile
#done

