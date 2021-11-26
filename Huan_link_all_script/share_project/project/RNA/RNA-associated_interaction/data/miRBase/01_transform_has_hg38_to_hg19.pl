# 提取hsa_hg38_to_hg19.gff3中的hg19,得文件01_hsa_hg19.gff3; 得转换hg19失败的hg38文件01_hsa_hg38_to_hg19_fail.gff3, gff3为1 based文件
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;


my $f1 = "./hsa_hg38_to_hg19.gff3";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n";
my $fo1 = "01_hsa_hg19.gff3"; #得转换hg19成功的文件
open my $O1, '>', $fo1 or die "$0 : failed to open output file  '$fo1' : $!\n";
my $fo2 = "01_hsa_hg38_to_hg19_fail.gff3"; #得转换hg19失败的hg38文件
open my $O2, '>', $fo2 or die "$0 : failed to open output file  '$fo2' : $!\n";
my $header = "Chr\t\tmiRNA_type\tPos1\tPos2\t\tStrand\t\tInfo";
print $O1 "$header\n";

while(<$I1>)
{
    chomp;
    if ($_=~/\->/){
        my @t = split/\->/;
        my $hg19 = $t[1];
        $hg19 =~ s/^\s+//g;
        print $O1 "$hg19\n";
    }
    else{
        print $O2 "$_\n";
    }
}
