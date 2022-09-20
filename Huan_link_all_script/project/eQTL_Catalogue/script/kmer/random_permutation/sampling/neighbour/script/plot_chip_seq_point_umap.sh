computeMatrix reference-point \
            -S $marker_bw  \
            -R ${cluster_inputdir}/cluster_1.bed \
                ${cluster_inputdir}/cluster_2.bed \
                ${cluster_inputdir}/cluster_3.bed \
                ${cluster_inputdir}/cluster_4.bed \
                ${cluster_inputdir}/cluster_5.bed \
                ${cluster_inputdir}/cluster_6.bed \
                ${cluster_inputdir}/cluster_7.bed \
                ${cluster_inputdir}/cluster_8.bed \
                ${cluster_inputdir}/cluster_9.bed \
                ${cluster_inputdir}/cluster_10.bed \
            --samplesLabel ${marker}\
            --referencePoint center \
            --binSize 50\
            -b 2500 -a 2500 \
            -out ${output_dir}/${marker}_${level}.gz
            # --skipZeros \


plotProfile -m ${output_dir}/${marker}_${level}.gz \
              -out ${output_dir}/${marker}_${level}.pdf\
              --colors "#1E77B4" "#FF7F0E" "#2CA02C" "#C22324" "#9567BD" "#8C554B" "#E277C1" "#7F7F7F" "#BCBC22" "#15BECE"\
              --plotHeight 10.5 \
              --plotWidth 10 \
              --refPointLabel 'center'  --regionsLabel "C1" "C2" "C3" "C4" "C5" "C6" "C7" "C8" "C9" "C10"\
              --yAxisLabel "" 

