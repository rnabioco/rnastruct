
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $DIR

echo  "
track rnastructurome
shortLabel rnastructurome
longLabel rnastucturome tracks
superTrack on show
"


METRICS=(
ed
fmfe
mfe
pvalue
zscore)


for metric in ${METRICS[@]}
do
    echo "    track RNA"$metric
    echo "    parent rnastructurome"
    echo "    shortName RNA"$metric
    echo "    longName rnastructurome $metric"
    echo "    visibility full" 
    echo "    maxHeightPixels 100:30:8"
    echo "    compositeTrack on"
    echo "    autoScale on"
    echo "    windowingFunction mean"
    echo "    type bigWig -100 100"  
    echo "    "  
    echo "    track rna"$metric"Pos" 
    echo "    bigDataUrl \
    http://amc-sandbox.ucdenver.edu/User33/bentley/hub/hg19/rnastructurome/"$metric"/all_chr_forward_"$metric"_hg19.bw"
    echo "    parent RNA"$metric 
    echo "    shortLabel rnastructurome"$metric"Pos"
    echo "    longLabel rnastructurome $metric positive" 
    echo "    autoScale on"
    echo "    type bigWig -100 100"  
    echo "    color 0,0,255"
    echo "    "
    echo "    track rna"$metric"Neg" 
    echo "    bigDataUrl \
    http://amc-sandbox.ucdenver.edu/User33/bentley/hub/hg19/rnastructurome/"$metric"/all_chr_reverse_"$metric"_hg19.bw"
    echo "    parent RNA"$metric
    echo "    shortLabel rnastructurome"$metric"Neg"
    echo "    longLabel rnastructurome $metric negative"
    echo "    type bigWig -100 100"  
    echo "    autoScale on"
    echo "    color 240,59,32" 
    echo ""
done
