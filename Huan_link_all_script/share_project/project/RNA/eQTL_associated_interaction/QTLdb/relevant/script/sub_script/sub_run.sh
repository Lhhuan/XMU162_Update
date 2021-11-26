#
ssh compute-0-4
cd /state/partition1/
mkdir sub_pre_file 
cd sub_pre_file
zless "/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/annotation_data/conservation/phastCons100way/merge_phastCons100way_sorted.bed.gz" | split -l 300000
cd ..
perl 01_bedtools_intersect.pl #将./sub_pre_file 中的文件与/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/output/ALL_eQTL/QTLbase_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_cis_trans_eQTL_normalized_sorted.bed.gz
#进行bedtools intersect，得文件存在./sub_intersect_file/下面
cat ./sub_intersect_file/* | gzip > QTLbase_all_eQTL_phastCons100way.bed.gz