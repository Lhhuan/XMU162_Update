perl 01_merge_all_eQTL_association_interaction.pl
echo -e "finish_01_merge_all_eQTL_association_interaction\n"
perl 02_add_pos_info_for_rsid.pl
echo -e "finish_02_add_pos_info_for_rsid\n"
perl 03_add_rsid.pl
echo -e "finish_03_add_rsid\n"
perl 04_add_ENSG_Entrez.pl
echo -e "finish_04_add_ENSG_Entrez\n"
perl 05_add_gene_position.pl
echo -e "finish_05_add_gene_position\n"
zless ../output/05_add_gene_position.txt.gz | cut -f2,3,21,22,24,25,26 >../output/05_used_clump.txt
# perl 06_plink_clump.pl
# echo -e "finish_06_plink_clump\n"