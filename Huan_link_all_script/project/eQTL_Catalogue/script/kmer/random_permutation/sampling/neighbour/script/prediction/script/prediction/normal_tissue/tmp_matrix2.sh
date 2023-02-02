computeMatrix reference-point \
    -S "/share/data0/QTLbase/huan/eQTL_Catalogue/prediction/normal_tissue/skin_suprapubic/marker/Human_FACTOR_merge_mean_signalvalue.bw" \
    -R "../predicted_region/predicted_regions_win5000_large_than6.bed"\
    --samplesLabel "Human_FACTOR"\
    --referencePoint center \
    --binSize 50\
    -b 2500 -a 2500 \
    --missingDataAsZero\
    -out ./tmp_out2/skin_suprapubic/Human_FACTOR_mean.bed.gz


computeMatrix reference-point \
    -S "/share/data0/QTLbase/huan/eQTL_Catalogue/prediction/normal_tissue/blood/marker/Human_FACTOR_merge_mean_signalvalue.bw" \
    -R "../predicted_region/predicted_regions_win5000_large_than6.bed"\
    --samplesLabel "Human_FACTOR"\
    --referencePoint center \
    --binSize 50\
    -b 2500 -a 2500 \
    --missingDataAsZero\
    -out ./tmp_out2/blood/Human_FACTOR_mean.bed.gz