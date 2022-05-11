#调整"/home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge_egene_0.05.bed.gz"格式得 "../output/01_all_hotspot.csv.gz"
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;


my $f1 = "../cell_line_info.txt";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n";
my $f2 = "../enhancer.txt";
open my $I2, '<', $f2 or die "$0 : failed to open input file '$f2' : $!\n";
# open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
my $fo1 = "../non_oncology_cell_line_list.txt"; #
# open my $O1, "| gzip >$fo1" or die $!;
open my $O1, '>', $fo1 or die "$0 : failed to open output file  '$fo1' : $!\n";
my $fo2 = "../enhancer_non_oncology_cell_line_list.txt"; #
open my $O2, '>', $fo2 or die "$0 : failed to open output file  '$fo2' : $!\n";
# print $O1 "Source node name,Target node name,Source node type,Target node type,Relationship type\n";

my (%hash1,%hash2);
while(<$I1>)
{
    chomp;
    my @f= split/\t/;
    for (my $i=0;$i<2;$i++){ #把所有未定义的空格等都替换成NONE
        unless(defined $f[$i]){
            $f[$i] = "NONE";
        }
    }
    unless(/^CELL_LINE/){
        my $cell= $f[0];
        my $cancer=$f[1];
        unless($cancer =~/YES/){
            print $O1 "$cell\n";
            $hash1{$cell}=1;
        }
    }
}


while(<$I2>)
{
    chomp;
    my @f= split/\t/;
    unless(/^CELL_LINE/){
        my $cell= $f[0];
        if(exists $hash1{$cell}){
            my $link= "http://www.enhanceratlas.org/data/download/enhancer/hs/${cell}.bed";
            system "wget -c $link";
            print "$cell\n";
            print $O2 "$cell\n";
        }
    }
}

