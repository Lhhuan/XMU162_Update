wget -c --no-check-certificate https://coxpresdb.jp/download/Hsa-r.c5-0/exp/Hsa-r.c5-0.expression.combat.txt.zip
wget -c --no-check-certificate https://coxpresdb.jp/download/Hsa-r.c5-0/coex/Hsa-r.v21-11.G16808-S60463.combat_pca.subagging.ls.d.zip #RNA-seq,

mv ../Hsa-r.v21-11.G16808-S60463.combat_pca.subagging.ls.d.zip ./G16808_S60463

unzip Hsa-r.v21-11.G16808-S60463.combat_pca.subagging.ls.d.zip

unzip Hsa-r.c5-0.expression.combat.txt.zip


wget -c --no-check-certificate https://coxpresdb.jp/download/Hsa-u.c3-0/coex/Hsa-u.v21-11.G16808-S85825.combat_pca.subagging.ls.d.zip #---RNA-seq and microarray 与webserer 展示一致，websever top300
mv Hsa-u.v21-11.G16808-S85825.combat_pca.subagging.ls.d.zip ./G16808_S85825 
cd G16808_S85825  
unzip Hsa-u.v21-11.G16808-S85825.combat_pca.subagging.ls.d.zip
mv Hsa-u.v21-11.G16808-S85825.combat_pca.subagging.ls.d.zip ../
cd ..
perl 01_get_gene_list.pl #获取./G16808_S85825 下的文件名，得01_G16808_S85825_get_gene_file_list.txt.gz
Rscript 02_transform_entrez_to_ENSG.R #02_G16808_S85825_entrezgene_symbol_ensembl.gene.txt
perl 03_unique_symbol_ENSG.pl ##unique "02_G16808_S85825_entrezgene_symbol_ensembl.gene.txt" 的ENSG,得03_entrezgene_symbol_ensembl.txt.gz
perl 04_transfrom_entrezgene_ensembl.pl #将 03_entrezgene_symbol_ensembl.txt.gz 中包含文件的top 300的entrezgene转成ensembl得./ENSG_G16808_S85825/${file_entrezgene}_${file_ensembl}