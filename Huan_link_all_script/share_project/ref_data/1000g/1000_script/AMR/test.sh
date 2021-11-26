#/bin/bash



for i in $( seq 22);
do
{
    bcftools view -S EUR_sample_list.txt /home/chaoqun/phase3/ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz -o 1kg.phase3.v5.shapeit2.eur.hg19.chr${i}.vcf.gz -O z
    sleep 5
    echo "bcftools view -S EUR_sample_list.txt /home/chaoqun/phase3/ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz -o 1kg.phase3.v5.shapeit2.eur.hg19.chr${i}.vcf.gz -O z"
} &
done
wait




#  for i in $( seq 22);