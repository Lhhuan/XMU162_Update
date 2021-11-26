for i in $( seq 22);
do
	echo $i;
    wget -c ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/phyloP100way/hg19.100way.phyloP100way/chr$i.phyloP100way.wigFix.gz
done