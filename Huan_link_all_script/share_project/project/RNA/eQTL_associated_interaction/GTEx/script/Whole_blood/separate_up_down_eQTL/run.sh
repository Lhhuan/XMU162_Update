perl 01_separate_up_down_eQTL_and_Completion_snp_for_xQTL_by_1kg.pl # # 用"/share/data0/1kg_phase3_v5_hg19/EUR/1kg.phase3.v5.shapeit2.eur.hg19.all.SNPs.vcf.gz" 补全"${dir}/${tissue}${suffix}"; 得up 补全文件${output_dir}/up_${tissue}_cis_eQTL_1kg_Completion.txt.gz #down 补全文件 ${output_dir}/down_${tissue}_cis_eQTL_1kg_Completion.txt.gz，得original 文件${output_dir}/up_${tissue}_cis_eQTL_original.txt.gz，${output_dir}/down_${tissue}_cis_eQTL_original.txt.gz
Rscript 02_up_NHP_big_par.R
Rscript 02_down_NHP_big_par.R
perl 03_filter_hotspot_interval18.pl #
perl 04_annotation_chromatin_states.pl ##interval_18 时，对"../../../output/${tissue}/separate_up_down_cis_eQTL/Cis_eQTL/${d_type}/hotspot/interval_18/${tissue}_segment_${group}_cutoff_${cutoff}.bed.gz"用annotation_chromatin_states_interval18.sh进行annotation,得$output_dir/$factor_$input_file_base_name
perl 05_filter_annotation_state.pl # 