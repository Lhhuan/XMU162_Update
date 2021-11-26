#统计../output/01_all_kinds_QTL.txt 中每特定距离中snp的数量，得../output/count_per_500bp_snp_number.txt, ../output/count_per_750bp_snp_number.txt, ../output/count_per_1000bp_snp_number.txt
#../output/count_per_500bp_snp_number_region.txt, ../output/count_per_750bp_snp_number_region.txt, ../output/count_per_1000bp_snp_number_region.txt

#!/usr/bin/perl
use warnings;
use strict; 
use utf8;


my @array = ("a","b","c","g");
my $oo = join("\t",@array[0..1]);
print "$oo\n";