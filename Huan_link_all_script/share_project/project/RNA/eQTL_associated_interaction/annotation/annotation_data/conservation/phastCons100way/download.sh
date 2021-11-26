for i in $( seq 22);
do
	echo $i;
    wget -c ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/phastCons100way/hg19.100way.phastCons/chr$i.phastCons100way.wigFix.gz
done