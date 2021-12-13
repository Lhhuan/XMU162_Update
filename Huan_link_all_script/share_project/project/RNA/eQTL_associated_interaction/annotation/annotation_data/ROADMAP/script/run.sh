perl 01_extract_roadmap_id.pl #download getx related marks to "/share/data0/GTEx/annotation/ROADMAP/sample/${roadmap_ID}
perl 02_merge_marker.pl ##将不同sample的相同mark进行合并，得 "/share/data0/GTEx/annotation/ROADMAP/sample/merge/${marker}_sorted_merge.bed.gz"

bedtools nuc -fi /public/reference/genome/hg38/hg38.fa -bed 200K.genome.3col >hg38.gcstat.txt