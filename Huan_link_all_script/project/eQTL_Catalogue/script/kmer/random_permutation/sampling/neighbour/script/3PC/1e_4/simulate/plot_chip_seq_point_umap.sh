computeMatrix reference-point \
            -S $marker_bw  \
            -R ${cluster_inputdir}/01_random_c1.bed \
                ${cluster_inputdir}/01_random_c2.bed \
                ${cluster_inputdir}/01_random_c3.bed \
            --samplesLabel ${marker}\
            --referencePoint center \
            --binSize 50\
            -b 2500 -a 2500 \
            -out ${output_dir}/${marker}_${level}.gz
            # --skipZeros \


plotProfile -m ${output_dir}/${marker}_${level}.gz \
              -out ${output_dir}/${marker}_${level}.pdf\
              --colors "#1E77B4" "#FF7F0E" "#2CA02C" \
              --plotHeight 10.5 \
              --plotWidth 10 \
              --refPointLabel 'center'  --regionsLabel "C1" "C2" "C3"\
              --yAxisLabel "" 

