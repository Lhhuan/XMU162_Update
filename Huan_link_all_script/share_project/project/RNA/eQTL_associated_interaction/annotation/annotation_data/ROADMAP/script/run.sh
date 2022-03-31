perl 01_extract_roadmap_id.pl #download getx related marks to "/share/data0/GTEx/annotation/ROADMAP/sample/${roadmap_ID}
perl 02_merge_marker.pl ##将不同sample的相同mark进行合并，得 "/share/data0/GTEx/annotation/ROADMAP/sample/merge/${marker}_sorted_merge.bed.gz"
mv /share/data0/GTEx/annotation/ROADMAP/sample/merge/ /share/data0/GTEx/annotation/ROADMAP/sample/GTEx_merge/
#-----------------以上为GTEx包括的biosample

#-----------eQTL_Catalogue_add_sample
perl 03_download_eQTL_Catalogue_add_sample.pl  #利用new_add_roadmap.txt 得eQTL_Catalogue_add_sample_id.txt  download eQTL_Catalogue_add_biosample to "/share/data0/GTEx/annotation/ROADMAP/sample/${roadmap_ID}
cat eQTL_Catalogue_add_sample_id.txt unique_roadmap_id.txt |sort -u >all_eQTL_Catalogue_sample_id_unique.txt
perl 04_merge_marker.pl #利用all_eQTL_Catalogue_sample_id_unique.txt将不同sample的相同mark进行合并，得 "/share/data0/GTEx/annotation/ROADMAP/sample/merge/${marker}_sorted_merge.bed.gz"
perl 05_tranform_hg19_to_hg38.pl  #将"/share/data0/GTEx/annotation/ROADMAP/sample/merge/${marker}_sorted_merge.bed.gz" 从 hg19转换到hg38,得"/share/data0/GTEx/annotation/ROADMAP/sample/merge/hg38/${marker}_sorted_merge.bed.gz"