perl 01_normal_miRNA_target.pl # 将../raw_data/miRTarBase_MTI.txt normal成 ../normalized/01_miRTarBase_normalized.txt
cat ../normalized/01_miRTarBase_normalized.txt | cut -f6,15| sort -u | wc -l #380640 -1
