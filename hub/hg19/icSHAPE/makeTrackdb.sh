#set directory to directory where bash script is located
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $DIR

echo "
track icShape 
superTrack on show
shortLabel icShape 
longLabel icShape Reactivities 

    track icShapeData
    parent icShape 
    compositeTrack on show
    shortLabel icShapeReactivities 
    longLabel icShape HEK293T Reactivites 
    type bigWig
    autoScale group 
    windowingFunction mean
    visibility full 
    maxHeightPixels 100:30:8
    subGroup1 source Source vivo=In_Vivo vitro=In_Vitro
    subGroup2 region Region cy=Cytoplasm np=Nucleoplasm ch=Chromatin
    subGroup3 strand Strand plus=Positive minus=Negative
    dimensions dimX=source dimY=region dimA=strand
    dimensionXchecked In_Vivo
    sortOrder region=+ source=-
"

color_pos="255,0,0"
color_neg="67,205,128"

for file in *.bw
  do 
      source_id=$(echo $file | cut -d "_" -f 3)
      region_id=$(echo $file | cut -d "_" -f 2)
      strand_id=$(echo $file | cut -d "_" -f 4 | sed 's/.bw//')

      trackname=$(basename $file .bw) 
      
      if [[ $strand_id == "minus" ]]
      then
          color=$color_neg
      else
          color=$color_pos
      fi
      
      if [[ $source_id == "vivo" && $region_id == "ch" ]]
      then
          vis="full"
      else
          vis="hide"
      fi

      echo "   track $trackname"
      echo "   parent icShapeData"
      echo "   shortLabel $trackname"
      echo "   bigDataUrl \
      http://amc-sandbox.ucdenver.edu/User33/bentley/hub/hg19/icSHAPE/$file"
      echo "   longLabel  icShape ($region_id $source_id $strand_id)"  
      echo "   subGroups source=$source_id region=$region_id strand=$strand_id"
      echo "   maxHeightPixels 30:30:10"
      echo "   visibility $vis"
      echo "   color $color"
      echo "   type bigWig"
      echo ""
 
done
