#set directory to directory where bash script is located
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $DIR

echo "
track DMS_Reactivity
superTrack on show
shortLabel DMS_reactivity 
longLabel Normalized DMS reactivity 

    track dms_bigwigs
    parent DMS_Reactivity 
    compositeTrack on show
    shortLabel DMS_reactivity_bigwigs 
    longLabel DMS reactivitity scores
    type bigWig
    autoScale on
    windowingFunction mean
    visibility hide 
    maxHeightPixels 100:30:8
    subGroup1 polII PolII WT=Wild_type C4=Slow_Mutant 
    subGroup2 norm Normalization_method norm=Global \
      bgsub=Background_subtracted raw=unnormalized 
    subGroup3 strand Strand pos=Positive neg=Negative
    dimensions dimX=polII dimY=norm dimA=strand
    dimensionXchecked Wild_type
    sortOrder strand=- polII=-
"
##

states=("WT" "C4") 
colors=("255,0,0" "67,205,128" "155,48,255" "25,25,112" "0,0,255" "255,165,0")
#prefix_for_state=("A_" "B_" "C_" "D_" "E_" "F_")
#set counter to grab color by index
color_idx=0

for state in ${states[@]}
do 
   
    #set counters to follow replicates
    which_rep=0
    med_rep=1
    hypo_rep=1
    fore_rep=1
    #get color value
    color=${colors[$color_idx]}
    
    #get prefix
    prefix=${prefix_for_state[$color_idx]}
    
    color_idx=$((color_idx + 1))
    
    for file in *$state*.bw
    do 
        raw=$(echo $file | grep -c "raw")
        bgsub=$(echo $file | grep -c "bgsub")
        norm=$(echo $file | grep -c "norm")
        pos_neg=$(echo $file | grep -c "pos")

        if [ $raw -eq 1 ]
        then
            method="raw"

        elif [ $bgsub -eq 1 ]
        then
            method="bgsub"
        
        elif [ $norm -eq 1 ]
        then    
            method="norm"
        
        else 
            echo "Filename not found to be medulla hypothamalus or brain the rest"  >&2
            exit 1
        
        fi


        if [ $pos_neg -eq 1 ]
        then 
            strand="pos"
        else
            strand="neg"
        fi
        
        trackname=$(basename $file .bw) 

        echo "track $trackname"
        echo "parent dms_bigwigs"
        echo "shortLabel $trackname"
        echo "bigDataUrl http://amc-sandbox.ucdenver.edu/User33/bentley/hub/hg19/dms/$file"
        echo "longLabel  $trackname"  
        echo "subGroups polII=$state norm=$method strand=$strand"
        echo "maxHeightPixels 30:30:10"
        echo "color $color"
        echo "type bigWig"
        echo ""
    
   done
done
