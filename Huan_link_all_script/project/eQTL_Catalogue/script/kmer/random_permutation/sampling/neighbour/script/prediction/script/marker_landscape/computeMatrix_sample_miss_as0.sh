computeMatrix reference-point \
            -S $marker_bw  \
            -R "./output/1_resemble_Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed" \
            --samplesLabel ${marker}\
            --referencePoint center \
            --binSize 50\
            -b 2500 -a 2500 \
            -p max\
            --missingDataAsZero\
            -out ${output_dir}/${marker}_${level}_negative.bed.gz
            # --skipZeros \

computeMatrix reference-point \
            -S $marker_bw  \
            -R "./output/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed" \
            --samplesLabel ${marker}\
            --referencePoint center \
            --binSize 50\
            -b 2500 -a 2500 \
            -p max\
            --missingDataAsZero\
            -out ${output_dir}/${marker}_${level}_positive.bed.gz
# plotProfile -m ${output_dir}/${marker}_${level}.gz \
#               -out ${output_dir}/${marker}_${level}.pdf\
#               --colors "#1E77B4" "#FF7F0E" \
#               --plotHeight 6.5 \
#               --plotWidth 6 \
#               --refPointLabel 'center'  --regionsLabel "C1" "C2" \
#               --yAxisLabel "" 

