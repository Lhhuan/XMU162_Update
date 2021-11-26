
bedtools makewindows -g ../output/find_max_emplambda.txt -w 1000 | gzip > ../output/interval_1kb.txt.gz
bedtools intersect -a ../output/interval_1kb.txt.gz  -b ../output/adjust_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_all_xQTL.bed.gz -wb -wa| gzip >../output/all_NHPoisson_emplambda_interval_1000_cutoff_7.3_1000_win.bed.gz

bedtools makewindows -g ../output/find_max_emplambda.txt -w 5000 | gzip > ../output/interval_5kb.txt.gz
bedtools intersect -a ../output/interval_5kb.txt.gz  -b ../output/adjust_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_all_xQTL.bed.gz -wb -wa| gzip >../output/all_NHPoisson_emplambda_interval_1000_cutoff_7.3_win_5kb_bed.gz

bedtools makewindows -g ../output/find_max_emplambda.txt -w 10000 | gzip > ../output/interval_10kb.txt.gz
bedtools intersect -a ../output/interval_10kb.txt.gz  -b ../output/adjust_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_all_xQTL.bed.gz -wb -wa| gzip >../output/all_NHPoisson_emplambda_interval_1000_cutoff_7.3_win_10kb_bed.gz

bedtools makewindows -g ../output/find_max_emplambda.txt -w 20000 | gzip > ../output/interval_20kb.txt.gz
bedtools intersect -a ../output/interval_20kb.txt.gz  -b ../output/adjust_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_all_xQTL.bed.gz -wb -wa| gzip >../output/all_NHPoisson_emplambda_interval_1000_cutoff_7.3_win_20kb_bed.gz

bedtools makewindows -g ../output/find_max_emplambda.txt -w 50000 | gzip > ../output/interval_50kb.txt.gz
bedtools intersect -a ../output/interval_50kb.txt.gz  -b ../output/adjust_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_all_xQTL.bed.gz -wb -wa| gzip >../output/all_NHPoisson_emplambda_interval_1000_cutoff_7.3_win_50kb_bed.gz

bedtools makewindows -g ../output/find_max_emplambda.txt -w 100000 | gzip > ../output/interval_100kb.txt.gz
bedtools intersect -a ../output/interval_100kb.txt.gz  -b ../output/adjust_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_all_xQTL.bed.gz -wb -wa| gzip >../output/all_NHPoisson_emplambda_interval_1000_cutoff_7.3_win_100kb_bed.gz

bedtools makewindows -g ../output/find_max_emplambda.txt -w 500000 | gzip > ../output/interval_500kb.txt.gz
bedtools intersect -a ../output/interval_500kb.txt.gz  -b ../output/adjust_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_all_xQTL.bed.gz -wb -wa| gzip >../output/all_NHPoisson_emplambda_interval_1000_cutoff_7.3_win_500kb_bed.gz

bedtools makewindows -g ../output/find_max_emplambda.txt -w 600000 | gzip > ../output/interval_600kb.txt.gz
bedtools intersect -a ../output/interval_600kb.txt.gz  -b ../output/adjust_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_all_xQTL.bed.gz -wb -wa| gzip >../output/all_NHPoisson_emplambda_interval_1000_cutoff_7.3_win_600kb_bed.gz

bedtools makewindows -g ../output/find_max_emplambda.txt -w 700000 | gzip > ../output/interval_700kb.txt.gz
bedtools intersect -a ../output/interval_700kb.txt.gz  -b ../output/adjust_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_all_xQTL.bed.gz -wb -wa| gzip >../output/all_NHPoisson_emplambda_interval_1000_cutoff_7.3_win_700kb_bed.gz

bedtools makewindows -g ../output/find_max_emplambda.txt -w 800000 | gzip > ../output/interval_800kb.txt.gz
bedtools intersect -a ../output/interval_800kb.txt.gz  -b ../output/adjust_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_all_xQTL.bed.gz -wb -wa| gzip >../output/all_NHPoisson_emplambda_interval_1000_cutoff_7.3_win_800kb_bed.gz

bedtools makewindows -g ../output/find_max_emplambda.txt -w 1000000 | gzip > ../output/interval_1MB.txt.gz
bedtools intersect -a ../output/interval_1MB.txt.gz  -b ../output/adjust_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_all_xQTL.bed.gz -wb -wa| gzip >../output/all_NHPoisson_emplambda_interval_1000_cutoff_7.3_win_1MB_bed.gz