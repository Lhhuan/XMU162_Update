#将./sub_pre_file 中的文件与/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/output/ALL_eQTL/QTLbase_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_cis_trans_eQTL_normalized_sorted.bed.gz
#进行bedtools intersect，得文件存在./sub_intersect_file/下面

#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use File::Basename;
# use Env qw(PATH);

my %hash1;


use warnings;
use strict;

my $dir = "./sub_pre_file";
my @file;
my $filename;


@file = `find $dir -type f`;

# my @aa;
# push @aa,$file[0];
# push @aa,$file[1];

# foreach $filename (@aa) {
foreach $filename (@file) {
    $filename =~s/\s+//g;
    my $f1 = $filename;
    my $f2 = basename($f1);

    # $ENV{'input_file'}  = $Normalized_sort_file;
    my $command = "bedtools intersect -a /home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/output/ALL_eQTL/QTLbase_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_cis_trans_eQTL_normalized_sorted.bed.gz -b $f1 -wa -wb > ./sub_intersect_file/$f2.bed";
    system $command;
    print "$command\n";
}