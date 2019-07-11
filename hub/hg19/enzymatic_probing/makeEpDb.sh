#set directory to directory where bash script is located
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $DIR

echo "
track EnzymaticProbing 
superTrack on show
shortLabel Enzymatic Probing 
longLabel Structure Score (Rnase Accessibility)

    track ep_bigwigs
    parent EnzymaticProbing 
    compositeTrack on show
    shortLabel Enzymatic_probing_bigwigs 
    longLabel RNase accessibility bigwigs
    type bigWig
    autoScale on
    windowingFunction mean
    visibility show 
    maxHeightPixels 100:30:8
    subGroup1 polII PolII WT=Wild_type C4=Slow_Mutant 
    subGroup2 strand Strand pos=Positive neg=Negative
    dimensions dimX=polII dimA=strand
    dimensionXchecked Wild_type
    sortOrder strand=- polII=-
"
##

states=("WT" "C4") 
colors=("155,48,255" "25,25,112" "0,0,255" "255,165,0")

#set counter to grab color by index
color_idx=0

for state in ${states[@]}
do 
   
    #get color value
    color=${colors[$color_idx]}
    
    #get prefix
    prefix=${prefix_for_state[$color_idx]}
    
    color_idx=$((color_idx + 1))
    
    for file in *$state*.bw
    do 
        pos_neg=$(echo $file | grep -c "pos")

        if [ $pos_neg -eq 1 ]
        then 
            strand="pos"
        else
            strand="neg"
        fi
        
        trackname=$(basename $file .bw) 

        echo "track $trackname"
        echo "parent ep_bigwigs" 
        echo "shortLabel $trackname"
        echo "bigDataUrl http://amc-sandbox.ucdenver.edu/User33/bentley/hub/hg19/enzymatic_probing/$file"
        echo "longLabel  $trackname"  
        echo "subGroups polII=$state strand=$strand"
        echo "maxHeightPixels 30:30:10"
        echo "color $color"
        echo "type bigWig"
        echo ""
    
   done
done
