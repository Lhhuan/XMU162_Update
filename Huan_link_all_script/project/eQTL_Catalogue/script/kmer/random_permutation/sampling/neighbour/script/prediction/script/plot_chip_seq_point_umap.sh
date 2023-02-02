computeMatrix reference-point \
            -S $marker_bw  \
            -R ${cluster_inputdir}/081_warm_region_predict_hotspot_true.bed \
                ${cluster_inputdir}/081_warm_region_predict_hotspot_false.bed \
            --samplesLabel ${marker}\
            --referencePoint center \
            --binSize 50\
            -b 2500 -a 2500 \
            -out ${output_dir}/${marker}_${level}.gz
            # --skipZeros \


plotProfile -m ${output_dir}/${marker}_${level}.gz \
              -out ${output_dir}/${marker}_${level}.pdf\
              --colors "#A593E0" "#F68657"\
              --plotHeight 6.5 \
              --plotWidth 6 \
              --refPointLabel 'center'  --regionsLabel "True" "False"\
              --yAxisLabel "" 

