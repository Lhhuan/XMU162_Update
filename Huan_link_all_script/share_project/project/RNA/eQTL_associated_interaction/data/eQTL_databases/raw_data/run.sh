#------------------ 
#MRCA eQTL database (Affymetrix expression array, n=405)
wget -c http://csg.sph.umich.edu/liang/eQTL/MRCA_5percentFDR.zip
unzip MRCA_5percentFDR.zip #得MRCA_5percentFDR.csv
#MRCE eQTL database (Illumina expression array,n=550)
wget -c http://csg.sph.umich.edu/liang/eQTL/MRCE_5percentFDR.csv.gz
#Meta-analysis results (_A indicates results from Affymetrix expression dataset and _I indicates results from Illumina expression dataset)
wget -c http://csg.sph.umich.edu/liang/eQTL/meta_5percentFDR.csv.gz
gzip MRCA_5percentFDR.csv