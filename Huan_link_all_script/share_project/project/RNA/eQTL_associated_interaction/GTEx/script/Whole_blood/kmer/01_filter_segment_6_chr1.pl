#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
use Parallel::ForkManager;

my $cutoff = 0.176;
my $group ="hotspot";
my $tissue = "Whole_Blood";


my $input_file = "../../../output/${tissue}/Cis_eQTL/${group}_cis_eQTL/interval_18/whole_blood_segment_${group}_cutoff_${cutoff}.bed.gz";
# my $output_file = "../../../output/ALL_${QTL}/cis_trans/interval_15/${group}/interval_${interval}_cutoff_7.3_${type}_${QTL}_segment_${group}.bed.gz";
my $input_file_base_name = basename($input_file);
my $output_dir = "../../../output/${tissue}/Cis_eQTL/${group}_cis_eQTL/interval_18_filter/6"; 
# mkdir $PMID;
# #------------
if(-e $output_dir){
    print "${output_dir}\texist\n";
}
else{
    system "mkdir -p $output_dir";
}
#------------

my $f1 =$input_file;
$input_file_base_name =~ s/whole_blood/Whole_Blood/g;
# open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
#----------------------
my $fo1 = "${output_dir}/${tissue}_segment_${group}_cutoff_${cutoff}_chr1.bed.gz";
# open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
open my $O1, "| gzip >$fo1" or die $!;
while(<$I1>)
{
    chomp;
    my @f= split/\t/;
    my $SNP_chr =$f[0];
    my $start =$f[1];
    my $end =$f[2];
    my $length = $end-$start;
    if($length >5 && $SNP_chr =~/\bchr1\b/ ){
        print $O1 "$_\n";
    }
}
print "$group\t$cutoff\n";
      