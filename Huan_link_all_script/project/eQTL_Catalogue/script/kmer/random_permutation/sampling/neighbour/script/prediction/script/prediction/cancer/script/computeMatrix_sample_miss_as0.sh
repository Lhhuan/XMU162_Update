computeMatrix reference-point \
            -S $marker_bw  \
            -R $input_file \
            --samplesLabel ${marker}\
            --referencePoint center \
            --binSize 50\
            -b 2500 -a 2500 \
            --missingDataAsZero\
            -out ${output_dir}/${marker}_${level}.bed.gz
            # --skipZeros \
