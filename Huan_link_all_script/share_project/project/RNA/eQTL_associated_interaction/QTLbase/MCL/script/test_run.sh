

# # perl 02_filter_hotspot_in_1kg_phase3_v5_eur.pl 
# # echo -e "02_filter_hotspot_in_1kg_phase3_v5_eur.pl\n"
# cd /state/partition1/huan/
# # plink --vcf /home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/MCL/output/02_hotspot_in_1kg_phase3_v5_eur.vcf.gz --recode --out 02_hotspot_in_1kg_phase3_v5_eur
# # echo -e "finish make file\n"
# # plink --file 02_hotspot_in_1kg_phase3_v5_eur --ld-window 1000 --ld-window-r2 0.2 --r2 --out 02_hotspot_in_1kg_phase3_v5_eur_1MB_100_0.2
# plink --file 02_hotspot_in_1kg_phase3_v5_eur --ld-window 10000 --ld-window-r2 0.2 --r2 --out 02_hotspot_in_1kg_phase3_v5_eur_1MB_10000_0.2
# echo -e "finish plink\n"
# gzip 02_hotspot_in_1kg_phase3_v5_eur_1MB_10000_0.2.ld
# echo -e "finish gzip\n"
# # gzip 02_hotspot_in_1kg_phase3_v5_eur_1MB_100_0.2.ld
# cd /state/partition1/huan/tmp/

# zcat ../02_hotspot_in_1kg_phase3_v5_eur_1MB_10000_0.2.ld.gz | split -l 10800000  -d -a 3 &&ls| grep x | xargs -n1 -i mv {} 02_hotspot_in_1kg_phase3_v5_eur_1MB_10000_0.2.ld_{}
# gzip *
# echo -e "split file\n"                                                 
# cd /home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/MCL/script/
perl 03_annotation_ld_band_par.pl 
echo -e "03_annotation_ld_band_par\n"
perl 031_cat_all_split_annotation.pl
echo -e "031_cat_all_split_annotation\n"