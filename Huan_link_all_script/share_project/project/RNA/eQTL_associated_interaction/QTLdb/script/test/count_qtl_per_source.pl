 #统计../output/all_qtl_clump_locus_r_square0.5.txt.gz中pop_tissue_qtl的数目，得../output/number_of_pop_tissue_qtl.txt,

#!/usr/bin/perl
use warnings;
use strict; 
use utf8;


my %hash1;

my $f1 = "../output/all_qtl_clump_locus_r_square0.5.txt.gz";
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
my $fo1 = "../output/number_of_pop_tissue_qtl.txt";
open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
my $header = "QTL_type\tTissue\tPopulation\tnumber";
print $O1 "$header\n";

while(<$I1>)
{
    chomp;
    unless(/PMID/){
        # print "$file\n";
        my @f =split/\t/;
        my $Population =$f[-1];
        my $Tissue =$f[-2];
        my $QTL_type =$f[-3];
        my $Locus_region =$f[-5];
        my $Tag_snp =$f[-6];
        my $k = "$QTL_type\t$Tissue\t$Population";
        my $v = "$Tag_snp\t$Locus_region";
        push @{$hash1{$k}},$v;
    }
}


foreach my $k (sort keys %hash1){
    my @vs =@{$hash1{$k}};
    my %hash2;
    @vs = grep { ++$hash2{$_} < 2 } @vs;
    my $number = @vs;
    print $O1 "$k\t$number\n";
}

close $O1;

system "cat ../output/number_of_pop_tissue_qtl.txt | sort -k4nr >../output/sorted_number_of_pop_tissue_qtl.txt";