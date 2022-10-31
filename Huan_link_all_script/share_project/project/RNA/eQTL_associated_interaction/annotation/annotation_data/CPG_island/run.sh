wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cpgIslandExt.txt.gz
zless cpgIslandExt.txt.gz |cut -f2-4 |sort -k1,1 -k2,2n |gzip >cpgIsland.bed.gz