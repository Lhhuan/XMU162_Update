computeMatrix reference-point \
            -S $marker_bw  \
            -R ${cluster_inputdir}/cluster_1.bed \
                ${cluster_inputdir}/cluster_2.bed \
                ${cluster_inputdir}/random_from_whole_genome.bed \
                ${cluster_inputdir}/random_from_cold_region.bed \
            --samplesLabel ${marker}\
            --referencePoint center \
            --binSize 50\
            -b 2500 -a 2500 \
            --missingDataAsZero\
            -out ${output_dir}/${marker}_${level}.gz
                        # -p max\
            # --skipZeros \


plotProfile -m ${output_dir}/${marker}_${level}.gz \
              -out ${output_dir}/${marker}_${level}.pdf\
              --colors "#1E77B4" "#FF7F0E" "#2CA02C" "#C22324"\
              --plotHeight 6.5 \
              --plotWidth 6 \
              --refPointLabel 'center'  --regionsLabel "C1" "C2" "Whole genome" "Cold region"\
              --yAxisLabel "" 

