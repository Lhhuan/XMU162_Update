plink --bfile /home/huanhuan/ref_data/1000g/1kg_phase3_v5_hg19/EUR/1kg.phase3.v5.shapeit2.eur.hg19.all.SNPs.uniq_posID --threads 20 --clump-r2 0.2 --clump-kb 500 --clump 1234.txt --out ./test/out_EUR
plink --bfile /home/huanhuan/ref_data/1000g/1kg_phase3_v5_hg19/AMR/1kg.phase3.v5.shapeit2.amr.hg19.all.SNPs.uniq.posID --threads 20 --clump-r2 0.2 --clump-kb 500 --clump 1234.txt --out ./test/out_AMR
#compute-0-3:
plink --bfile /state/partition1/huan/ref_data/1000g/1kg_phase3_v5_hg19/AFR/1kg.phase3.v5.shapeit2.amr.hg19.all.SNPs.uniq.posID --threads 20 --clump-r2 0.2 --clump-kb 500 --clump 1234.txt --out ./test/out_AFR
plink --bfile /state/partition1/huan/ref_data/1000g/1kg_phase3_v5_hg19/EAS/1kg.phase3.v5.shapeit2.eas.hg19.all.SNPs.uniq.posID --threads 20 --clump-r2 0.2 --clump-kb 500 --clump 1234.txt --out ./test/out_EAS
plink --bfile /state/partition1/huan/ref_data/1000g/1kg_phase3_v5_hg19/SAS/1kg.phase3.v5.shapeit2.sas.hg19.all.SNPs.uniq.posID --threads 20 --clump-r2 0.2 --clump-kb 500 --clump 1234.txt --out ./test/out_SAS




plink --bfile 1kg.phase3.v5.eur.hg19.SNPs.0.05 --clump-p1 1e-5 --clump-p2 0.01 --clump-r2 0.8 --clump-kb 250 --clump test.proxy.signal --clump-snp-field ID -clump-field PVAL --out test.proxy.signal



plink --bfile /share/data0/1kg_phase3_v5_hg19/EUR/1kg.phase3.v5.shapeit2.eur.hg19.all.SNPs.uniq_posID --show-tags demo/test.proxy --list-all --tag-r2 0.8

plink --bfile /share/data0/1kg_phase3_v5_hg19/EUR/1kg.phase3.v5.shapeit2.eur.hg19.all.SNPs.uniq_posID --clump-p1 1e-5 --clump-p2 0.01 --clump-r2 0.8 --clump-kb 250 --clump demo/test.proxy.signal --clump-snp-field ID -clump-field PVAL --clump-annotate --out demo/test.proxy.signal_anno

plink --bfile /share/data0/1kg_phase3_v5_hg19/EUR/1kg.phase3.v5.shapeit2.eur.hg19.all.SNPs.uniq_posID --clump-p1 1e-5 --clump-p2 0.01 --clump-r2 0.8 --clump-kb 250 --clump test.proxy.signal --clump-snp-field ID -clump-field PVAL --out test/test.proxy.signal
