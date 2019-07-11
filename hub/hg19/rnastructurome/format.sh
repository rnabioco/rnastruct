#BSUB -e err_%J.err 
#BSUB -o out_%J.out 
#BSUB -n 12 
#BSUB -J format[1-4] 
#BSUB -R "select[mem>10] rusage[mem=10] span[hosts=1]"


METRICS=(
ed
fmfe
mfe
pvalue
)

metric=${METRICS[$(($LSB_JOBINDEX - 1))]}

outdir=$metric

cd $outdir

loc="../hg38ToHg19.over.chain.gz"
chroms="../hg19.chrom.sizes"

bigWigToBedGraph \
    "all_chr_forward_"$metric".bw" stdout \
    | liftOver stdin $loc \
    "all_chr_forward_"$metric"_hg19.bg" tmp.txt

bigWigToBedGraph \
    "all_chr_reverse_"$metric".bw" stdout \
    | liftOver stdin $loc \
    "all_chr_reverse_"$metric"_hg19.bg" tmp.txt

sort --parallel=12 -S 8G -k1,1 -k2,2n \
    "all_chr_reverse_"$metric"_hg19.bg" \
    | bedtools merge -i - -d -1 -c 4 -o mean \
    > $metric"tmp_reverse.bg"

bedGraphToBigWig \
    $metric"tmp_reverse.bg" \
    $chroms \
    "all_chr_reverse_"$metric"_hg19.bw"


sort --parallel=12 -S 8G -k1,1 -k2,2n \
    "all_chr_forward_"$metric"_hg19.bg" \
    | bedtools merge -i - -d -1 -c 4 -o mean \
    > $metric"tmp_forward.bg"

bedGraphToBigWig \
    $metric"tmp_forward.bg" \
    $chroms \
    "all_chr_forward_"$metric"_hg19.bw"

rm $metric"tmp_forward.bg" $metric"tmp_reverse.bg"
