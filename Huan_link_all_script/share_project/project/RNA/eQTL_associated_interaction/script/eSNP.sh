#!/bin/bash
#PBS -N phase
#PBS -l nodes=1:compute-0-3:ppn=22,walltime=100000000000000:00:00
#PBS -o /home/chaoqun/phase3/err/phase.o
#PBS -e /home/chaoqun/phase3/err/phase.e
cd /home/chaoqun/phase3

#############################   cis
for i in $( seq 22);
do
	echo $i;
	for j in BRCA_N BRCA_P COAD KIRC LIHC LUAD OV PRAD STAD THCA UCEC;
	do
		awk -v a=chr$i 'BEGIN{print "SNP\tP"} {if($7==a)print $2"\t"$5}' /home/chaoqun/project/liwenzhi/MATRIX_EQTL/LNCRNA/eQTL/lncRNA_${j}_peer5_cis_anno.eqtl >/home/chaoqun/phase3/CLUMP/lncRNA_cis_${j}_${i}.eSNP
		plink --bfile $i --threads 20 --clump-r2 0.2 --clump-kb 500 --clump /home/chaoqun/phase3/CLUMP/lncRNA_cis_${j}_${i}.eSNP --out /home/chaoqun/phase3/CLUMP/lncRNA_cis_${j}_${i}
	done
done






########################    trans
for i in $( seq 22);
do
	echo $i;
	for j in BRCA_N BRCA_P COAD KIRC LIHC LUAD OV PRAD STAD THCA UCEC;
	do
		awk -v a=chr$i 'BEGIN{print "SNP\tP"} {if($7==a)print $2"\t"$5}' /home/chaoqun/project/liwenzhi/MATRIX_EQTL/LNCRNA/eQTL/lncRNA_${j}_peer5_trans_anno.eqtl >/home/chaoqun/phase3/CLUMP/lncRNA_trans_${j}_${i}.eSNP
		plink --bfile $i --threads 20 --clump-r2 0.2 --clump-kb 500 --clump /home/chaoqun/phase3/CLUMP/lncRNA_trans_${j}_${i}.eSNP --out /home/chaoqun/phase3/CLUMP/lncRNA_trans_${j}_${i}
	done
done


#for i in $( seq 22);
#do
#	plink --bfile ${i} \
#       --threads 20 \
#		--r2 \
#	--ld-snp-list ${i}.rs \
#	--ld-window-kb 1000 \
#	--ld-window 99999 \
#	--ld-window-r2 0.2 \
#	--out ${i}.LD
#done 
