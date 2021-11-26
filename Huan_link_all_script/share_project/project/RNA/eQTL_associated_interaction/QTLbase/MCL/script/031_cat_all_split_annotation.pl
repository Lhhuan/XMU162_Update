#将 "SNP_A\tSNP_B\tR2\tbind\n" 输出到 ../output/03_ALL_annotation_ld_band.txt.gz 并将 /state/partition1/huan/tmp_output/* >>../output/03_ALL_annotation_ld_band.txt.gz

#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;

my $fo1 = "../output/03_ALL_annotation_ld_band.txt.gz";
open my $O1, "| gzip >$fo1" or die $!;


print $O1  "SNP_A\tSNP_B\tR2\tbind\n";

close($O1);

my $command = "zcat /state/partition1/huan/tmp_output/* | gzip >>../output/03_ALL_annotation_ld_band.txt.gz";
print "$command\n";
system $command;
