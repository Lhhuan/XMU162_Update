computeMatrix reference-point \
            -S $marker_bw  \
            -R ${cluster_inputdir}/predicted_Class1.bed \
                ${cluster_inputdir}/predicted_Class2.bed \
                ${cluster_inputdir}/predicted_Class0.bed \
            --samplesLabel ${marker1}\
            --referencePoint center \
            --binSize 50\
            -b 2500 -a 2500 \
            --missingDataAsZero\
            -out ${output_dir}/${marker1}_${level}.gz
            # --skipZeros \
            # -p max\q

plotProfile -m ${output_dir}/${marker1}_${level}.gz \
              -out ${output_dir}/${marker1}_${level}.pdf\
              --colors "#1E77B4" "#FF7F0E" "#2CA02C" \
              --plotHeight 6.5 \
              --plotWidth 6 \
              --refPointLabel 'center'  --regionsLabel "C1" "C2" "C0"\
              --yAxisLabel "" 

