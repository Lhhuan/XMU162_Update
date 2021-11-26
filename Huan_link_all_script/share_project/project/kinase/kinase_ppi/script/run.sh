Rscript 01_transfrom_kinsae_gene_to_symbol.R # 将../data/kinase.txt 中的symbol 转为symbol（为了统一symbol），得../output/01_transfrom_kinsae_gene_to_symbol.txt
cat ../data/FIsInGene_122718_with_annotations.txt |awk 'NR>1'| cut -f1 | sort -u >../output/network_start.txt
cat ../data/FIsInGene_122718_with_annotations.txt |awk 'NR>1'| cut -f2 | sort -u >../output/network_end.txt
cat ../output/network_start.txt ../output/network_end.txt | sort -u >../output/all_network_gene.txt
Rscript 02_transfrom_network_gene_to_symbol.R #将 ../output/all_network_gene.txt 中的gene 转为symbol (为了得到统一symbol)，得../output/02_all_network_gene_symbol.txt
perl 03_add_no_transfrom_symbol_gene.pl #../output/all_network_gene.txt 中 没能转换成symbol的gene 用其原来的gene代替，
#得../output/03_final_all_network_gene_symbol.txt
perl 04_map_network_to_symbol.pl ##将../data/FIsInGene_122718_with_annotations.txt中的start和end用symbol表示（map到../output/03_final_all_network_gene_symbol.txt）
#得../output/04_network_symbol.txt
less ../output/04_network_symbol.txt | awk 'NR>1'| cut -f1,2,5 >../output/network_symbol_used_rwr.txt
perl 05_emerge_start_and_run_rwr_and_emerge_top.pl ##用../output/01_transfrom_kinsae_gene_to_symbol.txt 每个基因symbol作为成RWR的起点，并基于../output/network_symbol_used_rwr.txt
#run rwr，并把结果存在../output/RWR/rwr_result/文件夹下面
#取 rwr result的top 存在../output/RWR/rwr_result_top/文件夹下面
Rscript 06_transfrom_all_gene_to_symbol.R #将../data/mat_rank_gs_split_sig.txt gene转为symbol,得../output/06_all_gene_from_jintao_symbol.txt
perl 07_add_no_transfrom_symbol_gene.pl #将../data/mat_rank_gs_split_sig.txt中没有转成 symbol用原来的gene表示，
#得../output/07_final_all_gene_from_jintao_symbol.txt
perl 08_count_kinase_in_the_top_network.pl # 根据../output/07_final_all_gene_from_jintao_symbol.txt和../output/RWR/rwr_result_top/ 计算kinase in the top, 
#and kinase out of top, isn't no-kinase in the top and isn't no-kinase out of the top,得文件../output/08_fisher_need_data.txt
#得不在网络中的非kinase gene文件../output/08_nonkinase_not_in_the_network.txt 
#得在网络中的非kinase gene文件 ../output/08_nonkinase_in_the_network.txt
#得不在网络中的kinase gene文件../output/08_kinase_not_in_the_network.txt
#得在网络中的kinase gene文件 ../output/08_kinase_in_the_network.txt
#得symbol的$neighbour文件 ../output/08_start_neighbour.txt
#得 在不同start下的，neighbour得event文件 ../output/08_start_neighbour_events.txt
cat ../output/08_nonkinase_in_the_network.txt ../output/08_kinase_in_the_network.txt > ../output/08_all_gene_from_jintao_in_the_network.txt

Rscript 09_fisher_test.R #检测../output/08_fisher_need_data.txt的fisher 得../output/09_fisher_test_result.txt
perl 10_merge_fisher_result_and_events_neighbor.pl #将../output/09_fisher_test_result.txt和../output/08_start_neighbour_events.txt merge到一起
#得../output/10_start_neighbour_events_fisher.txt 
#得显著的文件../output/10_start_neighbour_events_fisher_significant.txt  
#得不显著的文件 ../output/10_start_neighbour_events_fisher_not_significant.txt  
#------------------------------------------------------------------------------#查看rwr结果的分布
perl 12_filter_kinase_from_network_and_random_select_start.pl #将../output/03_final_all_network_gene_symbol.txt 中的../output/08_kinase_in_the_network.txt
#过滤掉，并随机选择69个基因。得../output/12_random_start.txt
perl 13_random_start_rwr.pl #用得../output/12_random_start.txt的单个基因做为start，用../output/network_symbol_used_rwr.txt作为网络，走rwr，
#结果存在../RWR_random/rwr_result/文件夹下面,取rwr result的top 存在../output/RWR_random/rwr_result_top/文件夹下面
cat ../output/RWR/rwr_result/*.txt | cut -f2 >../output/all_kinase_rwr_result.txt
cat ../output/RWR_random/rwr_result/*.txt | cut -f2 >../output/all_random_rwr_result.txt

# perl 14_merge_all_kinase_and_random_rwr_result.pl # 把../output/RWR/rwr_result/*.


#---------------------------------------------------------------------用jintao所有的gene 走random walk并取top
perl 11_emerge_start_and_run_rwr_and_emerge_top.pl ###../output/08_all_gene_from_jintao_in_the_network.txt 每个基因symbol作为成RWR的起点，并基于../output/network_symbol_used_rwr.txt
#run rwr，并把结果存在../output/RWR_a_g_JT/rwr_result/文件夹下面
#取 rwr result的top 存在../output/RWR_a_g_JT/rwr_result_top/文件夹下面
perl 12_filter_kinase_from_network_and_random_select_start.pl ##将../output/03_final_all_network_gene_symbol.txt 中的../output/08_kinase_in_the_network.txt
#过滤掉，并随机选择69个基因。得../output/12_random_start.txt
perl 13_random_start_rwr.pl ##用得../output/12_random_start.txt的单个基因做为start，用../output/network_symbol_used_rwr.txt作为网络，走rwr，
#结果存在../RWR_random/rwr_result/文件夹下面,取rwr result的top 存在../output/RWR_random/rwr_result_top/文件夹下面
Rscript 14_density_polt_merge.R #







perl 06_count_kinase_in_the_top_network.pl #根据../output/01_transfrom_kinsae_gene_to_symbol.txt和../output/RWR/rwr_result_top/ 计算kinase in the top, 
#and kinase out of top, isn't no-kinase in the top and isn't no-kinase out of the top,得文件../output/06_fisher_need_data.txtss
#得不在网络中的kinase文件../output/06_kinase_not_in_the_network.txt