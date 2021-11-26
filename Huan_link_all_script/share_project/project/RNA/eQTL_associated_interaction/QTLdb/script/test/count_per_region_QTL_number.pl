#统计../output/all_qtl_clump_locus_r_square0.5.txt.gz 中每特定距离中QTL(tag snp的数量)，得../output/count_per_100kb_QTL_number.txt
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;


my $f1 = "../output/all_qtl_clump_locus_r_square0.5.txt.gz";
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
my $fo1 = "../output/count_per_100kb_QTL_number.txt";
open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";

my $header = "chr\tregion_start\tregion_end\tnumber";
print $O1 "$header\n";

my (%hash1, %hash2, %hash3);

while(<$I1>){
    chomp;
    unless(/SNP/){
        # print "$file\n";
        my @f =split/\t/;
        my $Population =$f[-1];
        my $Tissue =$f[-2];
        my $QTL_type =$f[-3];
        my $Locus_region =$f[-5];
        my $Tag_snp =$f[-6];
        if ($Population =~/EUR/ && $Tissue =~ /\bBlood\b/ && $QTL_type =~/eQTL/){
            my @t = split/\_/,$Tag_snp;
            my $chr = $t[0];
            my $pos = $t[1];
            my $time1 = int($pos/100000); #
            my $k1 = "$chr\t$time1";
            push @{$hash1{$k1}},$Tag_snp;
        }
    }
}
#-----------------------100kb
foreach my $k (sort keys %hash1){ 
    my @v = @{$hash1{$k}};
    my %hash4;
    @v = grep { ++$hash4{$_} < 2 } @v;
    my $number =@v;
    my @t =split/\t/,$k;
    my $chr =$t[0];
    my $time = $t[1];
    my $pos_int = $time *100000;
    my $pos_end = $pos_int + 100000;
    my $region = "${pos_int}_${pos_end}";
    print $O1 "$chr\t$pos_int\t$pos_end\t$number\n";
}
