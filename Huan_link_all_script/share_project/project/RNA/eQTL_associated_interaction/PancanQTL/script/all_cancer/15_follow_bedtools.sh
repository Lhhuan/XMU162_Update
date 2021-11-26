cd ../../output/cancer_total/specific/pure/THCA
gunzip cancer_THCA_Thyroid_specific.bed.gz
bedtools intersect -a cancer_THCA_Thyroid_specific.bed -b qtl_gene.bed.gz -wa -wb |gzip  >cancer_THCA_Thyroid_specific_gene.bed.gz
gzip cancer_THCA_Thyroid_specific.bed