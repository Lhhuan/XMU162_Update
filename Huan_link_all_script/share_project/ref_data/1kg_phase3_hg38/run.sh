wget -c http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/phased-manifest_July2021.tsv
perl 01_filter_and_download.pl
wget -c https://www.internationalgenome.org/data-portal/data-collection/30x-grch38/igsr-1000 genomes 30x on grch38.tsv.tsv #from website rename igsr_1000_genomes_30x_on_grch38.tsv
perl 02_filter_EUR_sample_name.pl


#perl 02_extract_EUR_samples_in_vcf_format.pl  中多线程和bash test.sh 多线程 提EUR局费时且不准确，so,手动多线程
nohup bcftools view -S EUR_sample_list.txt "20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr1.recalibrated_variants.vcf.gz" -o EUR_20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr1.recalibrated_variants.vcf.gz -O z &
nohup bcftools view -S EUR_sample_list.txt "20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr2.recalibrated_variants.vcf.gz" -o EUR_20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr2.recalibrated_variants.vcf.gz -O z &
nohup bcftools view -S EUR_sample_list.txt "20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr3.recalibrated_variants.vcf.gz" -o EUR_20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr3.recalibrated_variants.vcf.gz -O z &
nohup bcftools view -S EUR_sample_list.txt "20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr4.recalibrated_variants.vcf.gz" -o EUR_20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr4.recalibrated_variants.vcf.gz -O z &
nohup bcftools view -S EUR_sample_list.txt "20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr5.recalibrated_variants.vcf.gz" -o EUR_20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr5.recalibrated_variants.vcf.gz -O z &
nohup bcftools view -S EUR_sample_list.txt "20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr6.recalibrated_variants.vcf.gz" -o EUR_20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr6.recalibrated_variants.vcf.gz -O z &
nohup bcftools view -S EUR_sample_list.txt "20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr7.recalibrated_variants.vcf.gz" -o EUR_20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr7.recalibrated_variants.vcf.gz -O z &
nohup bcftools view -S EUR_sample_list.txt "20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr8.recalibrated_variants.vcf.gz" -o EUR_20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr8.recalibrated_variants.vcf.gz -O z &
nohup bcftools view -S EUR_sample_list.txt "20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr9.recalibrated_variants.vcf.gz" -o EUR_20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr9.recalibrated_variants.vcf.gz -O z &
nohup bcftools view -S EUR_sample_list.txt "20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr10.recalibrated_variants.vcf.gz" -o EUR_20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr10.recalibrated_variants.vcf.gz -O z &
nohup bcftools view -S EUR_sample_list.txt "20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr11.recalibrated_variants.vcf.gz" -o EUR_20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr11.recalibrated_variants.vcf.gz -O z &
nohup bcftools view -S EUR_sample_list.txt "20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr12.recalibrated_variants.vcf.gz" -o EUR_20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr12.recalibrated_variants.vcf.gz -O z &
nohup bcftools view -S EUR_sample_list.txt "20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr13.recalibrated_variants.vcf.gz" -o EUR_20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr13.recalibrated_variants.vcf.gz -O z &
nohup bcftools view -S EUR_sample_list.txt "20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr14.recalibrated_variants.vcf.gz" -o EUR_20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr14.recalibrated_variants.vcf.gz -O z &
nohup bcftools view -S EUR_sample_list.txt "20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr15.recalibrated_variants.vcf.gz" -o EUR_20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr15.recalibrated_variants.vcf.gz -O z &
nohup bcftools view -S EUR_sample_list.txt "20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr16.recalibrated_variants.vcf.gz" -o EUR_20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr16.recalibrated_variants.vcf.gz -O z &nohup bcftools view -S EUR_sample_list.txt "20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr17.recalibrated_variants.vcf.gz" -o EUR_20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr17.recalibrated_variants.vcf.gz -O z &
nohup bcftools view -S EUR_sample_list.txt "20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr18.recalibrated_variants.vcf.gz" -o EUR_20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr18.recalibrated_variants.vcf.gz -O z &
nohup bcftools view -S EUR_sample_list.txt "20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr19.recalibrated_variants.vcf.gz" -o EUR_20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr19.recalibrated_variants.vcf.gz -O z &
nohup bcftools view -S EUR_sample_list.txt "20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr20.recalibrated_variants.vcf.gz" -o EUR_20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr20.recalibrated_variants.vcf.gz -O z &
nohup bcftools view -S EUR_sample_list.txt "20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr21.recalibrated_variants.vcf.gz" -o EUR_20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr21.recalibrated_variants.vcf.gz -O z &
nohup bcftools view -S EUR_sample_list.txt "20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr22.recalibrated_variants.vcf.gz" -o EUR_20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr22.recalibrated_variants.vcf.gz -O z &