#统计../output/unique_by_qtl_tissue_pop_locus.txt.gz中不同区段长度的数目，得../output/count_of_different_distance_of_region.txt
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;

my $f1 = "../output/unique_by_qtl_tissue_pop_locus.txt.gz";
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件

my $fo1 = "../output/count_of_different_distance_of_region_QTL_specific.txt";
open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
my $fo2 = "../output/count_of_different_distance_of_region.txt";
open my $O2, '>', $fo2 or die "$0 : failed to open output file '$fo2' : $!\n";

print $O1 "QTL_type\tdistance\tnumber\n";
print $O2 "distance\tnumber\n";


my %hash1;
my %hash2;

# my @cutoffs = ("0", "20000","40000","60000","80000","100000","200000","300000","400000","500000");

while(<$I1>)
{
    chomp;
    unless(/^Locus_region/){
        my @f =split/\t/;
        my $Locus_region = $f[0];
        my $Locus_distance =$f[1];
        my $QTL_type = $f[2];
        my $Tissue =$f[3];
        my $Population =$f[4];
        my $v = $Locus_region;
        if ($Locus_distance >=0 && $Locus_distance<100000){
                my $k1_a = "$QTL_type\t0-100000";
                my $k2_a = "0-100000";
                push @{$hash1{$k1_a}},$v;
                push @{$hash2{$k2_a}},$v;
            if ($Locus_distance >=0 && $Locus_distance<20000){ #0-20000
                my $k1 = "$QTL_type\t0-20000";
                my $k2 = "0-20000";
                push @{$hash1{$k1}},$v;
                push @{$hash2{$k2}},$v;
            }
            elsif($Locus_distance >=20000 && $Locus_distance<40000){
                my $k1 = "$QTL_type\t20000-40000";
                my $k2 = "20000-40000";
                push @{$hash1{$k1}},$v;
                push @{$hash2{$k2}},$v;            
            }
            elsif($Locus_distance >=40000 && $Locus_distance<60000){
                my $k1 = "$QTL_type\t40000-60000";
                my $k2 = "40000-60000";
                push @{$hash1{$k1}},$v;
                push @{$hash2{$k2}},$v;            
            }
            elsif($Locus_distance >=60000 && $Locus_distance<80000){
                my $k1 = "$QTL_type\t60000-80000";
                my $k2 = "60000-80000";
                push @{$hash1{$k1}},$v;
                push @{$hash2{$k2}},$v;            
            }
            elsif($Locus_distance >=80000 && $Locus_distance<100000){
                my $k1 = "$QTL_type\t80000-100000";
                my $k2 = "80000-100000";
                push @{$hash1{$k1}},$v;
                push @{$hash2{$k2}},$v;            
            }
        }
        elsif($Locus_distance >=100000 && $Locus_distance<200000){
            my $k1 = "$QTL_type\t100000-200000";
            my $k2 = "100000-200000";
            push @{$hash1{$k1}},$v;
            push @{$hash2{$k2}},$v;            
        }
        elsif($Locus_distance >=200000 && $Locus_distance<300000){
            my $k1 = "$QTL_type\t200000-300000";
            my $k2 = "200000-300000";
            push @{$hash1{$k1}},$v;
            push @{$hash2{$k2}},$v;            
        }
        elsif($Locus_distance >=300000 && $Locus_distance<400000){
            my $k1 = "$QTL_type\t300000-400000";
            my $k2 = "300000-400000";
            push @{$hash1{$k1}},$v;
            push @{$hash2{$k2}},$v;            
        }
        else{
            my $k1 = "$QTL_type\t400000-500000";
            my $k2 = "400000-500000";
            push @{$hash1{$k1}},$v;
            push @{$hash2{$k2}},$v;  
        }
    }
}




foreach my $k(sort keys %hash1){
    my @v1s= @{$hash1{$k}};
    my %hash3;
    @v1s = grep { ++$hash3{$_} < 2 } @v1s;
    my $number = @v1s;
    print $O1 "$k\t$number\n";
}

foreach my $k(sort keys %hash2){
    my @v2s= @{$hash2{$k}};
    my %hash3;
    @v2s = grep { ++$hash3{$_} < 2 } @v2s;
    my $number = @v2s;
    print $O2 "$k\t$number\n";
}
