 #将04_all_cell_line_info.txt 中对应的normal cell line 从/home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/CTCF/data/中选出来，相比于05_filter_normal_cell_line.pl的输出多sample一列，得05_normal_cell_line_ctcf_sample.bed.gz,多输出signalvalue得05_normal_cell_line_ctcf_signalvalue.bed.gz


#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
use Parallel::ForkManager;

my $f1 = "04_all_cell_line_info.txt";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
my $f2 = "../data/download.sh";
open my $I2, '<', $f2 or die "$0 : failed to open input file '$f2' : $!\n"; 
# open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
#-------------------
my $fo1 = "./05_normal_cell_line_ctcf_sample.bed.gz";
open my $O1, "| gzip >$fo1" or die $!;
my $fo2 = "./05_normal_cell_line_ctcf_signalvalue.bed.gz";
open my $O2, "| gzip >$fo2" or die $!;
# open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";

# print $O2 "Cell_line\tdisease\n";

my %hash1;
my %hash2;

while(<$I1>){
    chomp;
    my @f = split/\t/;
    unless(/^Cell_line/){
        my $CELL_Line =$f[0];
        my $Disease =$f[1];
        if ($Disease =~/\bNO\b/){
            $hash1{$CELL_Line}=$Disease;
        }
    }
}
my %hash4;
for (my $i=1;$i<23;$i++){
    my $k = "chr${i}";
    $hash4{$k}=1;
}

while(<$I2>){
    chomp;
    my @f = split/\s+/;
    my $link = $f[2];
    my @t = split/\//,$link;
    # $link =~ s///g;
    # print  "$t[7]\n"; 
    my $org_file = $t[7];
    my $file_name = $org_file;
    $file_name =~ s/wgEncodeUwTfbs//g;   
    $file_name =~ s/Ctcf.*//g;   
    # wgEncodeUwTfbsA549CtcfStdPkRep2.narrowPeak.gz
    if(exists $hash1{$file_name}){
        my $f3 = "../data/$org_file";
        print  "$f3\n";
        # open my $I3, '<', $f3 or die "$0 : failed to open input file '$f3' : $!\n"; 
        open( my $I3 ,"gzip -dc $f3|") or die ("can not open input file '$f3' \n"); #读压缩文件
        while(<$I3>)
        {
            chomp;
            # unless(/^chr/){
            my @f = split/\s+/;
            my $chr = $f[0];
            my $start =$f[1];
            my $end =$f[2];
            my $signalvalue =$f[6];
            my $output1 = "$chr\t$start\t$end\t$file_name"; #cell line
            my $output2 = "$chr\t$start\t$end\t$signalvalue";
            if (exists $hash4{$chr}){
                print $O2 "$output2\n";
                unless(exists $hash2{$output1}){
                    $hash2{$output1}=1;
                    print $O1 "$output1\n";
                }
            }
        }
    }
    
    
}