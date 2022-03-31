bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%ALT'  "/share/data0/1kg_phase3_v5_hg19/EUR/1kg.phase3.v5.shapeit2.eur.hg19.all.vcf.gz" -O z -o 1kg.phase3.v5.shapeit2.eur.hg19.all.posID.vcf.gz
bcftools query -f '%ID\t%EUR_AF\n' 1kg.phase3.v5.shapeit2.eur.hg19.all.posID.vcf.gz |gzip >  1kg.phase3.v5.shapeit2.eur.hg19.all.maf.id.vcf.gz 
perl 05_filter_id_by_maf_greater0.pl

