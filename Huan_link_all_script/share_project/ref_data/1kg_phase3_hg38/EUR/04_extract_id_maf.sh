for((i=1;i<=22;i++));  
do   
    echo $i;  
    bcftools query -f '%ID\t%AF_EUR\n' EUR_2020-08-05_chr${i}.filtered.shapeit2-duohmm-phased.vcf.gz |gzip >  ../ID_MAF/EUR_2020-08-05_chr${i}.filtered.shapeit2-duohmm-phased.maf.id.txt.gz  
done  