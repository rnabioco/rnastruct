chroms="/beevol/home/riemondy/Projects/rnastruct/pipeline/norm_enzymatic_data/chroms.txt"

bedGraphToBigWig \
    C4_footprint_sorted_expressed_10.bedgraph \
    $chroms \
    C4_footprint_sorted_expressed_10.bw


bedGraphToBigWig \
    WT_footprint_sorted_expressed_10.bedgraph \
    $chroms \
    WT_footprint_sorted_expressed_10.bw


