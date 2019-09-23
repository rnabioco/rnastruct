#BSUB -n 1
#BSUB -J filter[1-4]
#BSUB -e err_%J.txt
#BSUB -e out_%J.txt
#BSUB -R "rusage[mem=10] select[mem>10] span[hosts=1]"
#BSUB -q rna


DEPTHS=(
10
15
50
100
)

depth=${DEPTHS[$(($LSB_JOBINDEX - 1))]}
chroms="/beevol/home/riemondy/Projects/shared_dbases/genomes/human/chr_appended_grch37/Homo_sapiens.GRCh37.dna.primary_assembly_chrappended.fa.fai"

# keep regions with 40% of nucleotides at coverage cutoff 
# only counting A and C
# slop by 10bp for folding

mkdir -p no_slop

bedtools map \
    -c 5 \
    -o count \
    -s \
    -a evofold.bed \
    -b "../../data/mismatches/july2019/depth_beds/"$depth"/WT/depth.bg" \
  | awk '$10 > 0 {print $0,$10/$7}' \
  | awk '$11>0.40' \
  > "no_slop/"$depth"_cutoff.txt"

