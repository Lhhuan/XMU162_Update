#统计../output/01_all_kinds_QTL.txt 中每特定距离中snp的数量，得../output/count_per_500bp_snp_number.txt, ../output/count_per_750bp_snp_number.txt, ../output/count_per_1000bp_snp_number.txt
#../output/count_per_500bp_snp_number_region.txt, ../output/count_per_750bp_snp_number_region.txt, ../output/count_per_1000bp_snp_number_region.txt

#!/usr/bin/perl
use warnings;
use strict; 
use utf8;


my (%hash1, %hash2, %hash3,%hash5);
my $f1 = "../data/QTLbase_download_data_sourceid.txt";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
my $f2 = "../output/01_all_kinds_QTL.txt";
open my $I2, '<', $f2 or die "$0 : failed to open input file '$f2' : $!\n"; 
my $fo1 = "../output/count_per_500bp_snp_number.txt";
open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
my $fo2 = "../output/count_per_750bp_snp_number.txt";
open my $O2, '>', $fo2 or die "$0 : failed to open output file '$fo2' : $!\n";
my $fo3 = "../output/count_per_1000bp_snp_number.txt";
open my $O3, '>', $fo3 or die "$0 : failed to open output file '$fo3' : $!\n";
my $fo4 = "../output/count_per_500bp_snp_number_region.txt";
open my $O4, '>', $fo4 or die "$0 : failed to open output file '$fo4' : $!\n";
my $fo5 = "../output/count_per_750bp_snp_number_region.txt";
open my $O5, '>', $fo5 or die "$0 : failed to open output file '$fo5' : $!\n";
my $fo6 = "../output/count_per_1000bp_snp_number_region.txt";
open my $O6, '>', $fo6 or die "$0 : failed to open output file '$fo6' : $!\n";

my $header = "chr\tregion_start\tregion_end\tnumber";
print $O1 "$header\n";
print $O2 "$header\n";
print $O3 "$header\n";
my $header2 = "chr\tregion\tnumber";
print $O4 "$header2\n";
print $O5 "$header2\n";
print $O6 "$header2\n";


while(<$I1>)
{
    chomp;
    unless(/PMID/){
        # print "$file\n";
        my @f =split/\t/;
        my $Sourceid = $f[1];
        my $Tissue =$f[3];
        my $Population = $f[4];
        my $v = "$Tissue\t$Population";
        $hash5{$Sourceid}=$v; #用于为$f2添加$Tissue和$Population
    }
}


while(<$I2>)
{
    chomp;
    unless(/SNP_chr/){
        # print "$file\n";
        my @f =split/\t/;
        my $chr = $f[0];
        my $pos = $f[1];
        my $Sourceid =$f[-2];
        my $QTL_type =$f[-1];
        my $snp = "${chr}_${pos}";
        if (exists $hash5{$Sourceid}){
            my $v = $hash5{$Sourceid};
            my @t = split/\t/,$v;
            my $Tissue =$t[0];
            my $Population =$t[1];
            if ($Population =~/EUR/ && $Tissue =~ /\bBlood\b/ && $QTL_type =~/eQTL/){
                my $time1 = int($pos/500); #
                # print "$time1\n";
                my $time2 = int($pos/750);
                my $time3 = int($pos/1000);
                my $k1 = "$chr\t$time1";
                my $k2 = "$chr\t$time2";
                my $k3 = "$chr\t$time3";
                push @{$hash1{$k1}},$snp;
                push @{$hash2{$k2}},$snp;
                push @{$hash3{$k3}},$snp;
            }
        }
    }
}
#-----------------------500bp
foreach my $k (sort keys %hash1){ 
    my @v = @{$hash1{$k}};
    my %hash4;
    @v = grep { ++$hash4{$_} < 2 } @v;
    my $number =@v;
    my @t =split/\t/,$k;
    my $chr =$t[0];
    my $time = $t[1];
    my $pos_int = $time *500;
    my $pos_end = $pos_int + 500;
    my $region = "${pos_int}_${pos_end}";
    print $O1 "$chr\t$pos_int\t$pos_end\t$number\n";
    print $O4 "$chr\t$region\t$number\n";
}

#------------------------750bp
foreach my $k (sort keys %hash2){ 
    my @v = @{$hash2{$k}};
    my %hash4;
    @v = grep { ++$hash4{$_} < 2 } @v;
    my $number =@v;
    my @t =split/\t/,$k;
    my $chr =$t[0];
    my $time = $t[1];
    my $pos_int = $time *750;
    my $pos_end = $pos_int + 750;
    my $region = "${pos_int}_${pos_end}";
    print $O2 "$chr\t$pos_int\t$pos_end\t$number\n";
    print $O5 "$chr\t$region\t$number\n";
}

#-------------------------------

foreach my $k (sort keys %hash3){ 
    my @v = @{$hash3{$k}};
    my %hash4;
    @v = grep { ++$hash4{$_} < 2 } @v;
    my $number =@v;
    my @t =split/\t/,$k;
    my $chr =$t[0];
    my $time = $t[1];
    my $pos_int = $time *1000;
    my $pos_end = $pos_int + 1000;
    my $region = "${pos_int}_${pos_end}";
    print $O3 "$chr\t$pos_int\t$pos_end\t$number\n";
    print $O6 "$chr\t$region\t$number\n";
}