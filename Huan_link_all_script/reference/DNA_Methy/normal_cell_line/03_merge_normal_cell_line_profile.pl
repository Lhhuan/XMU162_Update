#将02_all_cell_line_info.txt 中正常细胞的bed合并，得03_normal_cell_line_profile.bed

#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
use Parallel::ForkManager;

my $f1 = "02_all_cell_line_info.txt";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
my $f2 = "../01_need_file.txt";
open my $I2, '<', $f2 or die "$0 : failed to open input file '$f2' : $!\n"; 

#-------------------
my $fo1 = "./03_normal_cell_file.txt";
# open my $O1, "| gzip >$fo1" or die $!;
open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
my $fo3 = "./03_normal_cell.txt";
# open my $O1, "| gzip >$fo1" or die $!;
open my $O3, '>', $fo3 or die "$0 : failed to open output file '$fo3' : $!\n";
my $fo2 = "./03_normal_cell_line_profile.bed.gz";
open my $O2, "| gzip >$fo2" or die $!;
# open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
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

while(<$I2>)
{
    chomp;
    my @f = split/\s+/;
    my $file_name = $f[0];
    my $cell =$f[7];
    my $treatment=$f[8];
    my $file_type= $f[-3];
    my $size =$f[-1];
    $cell =~s/cell=|;//g;
    if(exists $hash1{$cell}){
        print $O1 "$_\n";
        print $O3 "$cell\n";
        my $f3 = "../raw_data/$file_name"; # 格式参考 https://genome.ucsc.edu/cgi-bin/hgTables?db=hg19&hgta_group=regulation&hgta_track=wgEncodeHaibMethylRrbs&hgta_table=wgEncodeHaibMethylRrbsBcskeletalmuscleh12817nBiochainSitesRep1&hgta_doSchema=describe+table+schema
        # open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
        open( my $I3 ,"gzip -dc $f3|") or die ("can not open input file '$f3' \n"); #读压缩文件
        while(<$I3>){
            chomp;
            my @f = split/\t/;
            unless(/^track/){
                my $output =join("\t",@f[0..2]);
                my $info=join(";",@f[3,4,9,10]);
                print $O2 "$output\t$info\n";
            }
        }

    }
    
}