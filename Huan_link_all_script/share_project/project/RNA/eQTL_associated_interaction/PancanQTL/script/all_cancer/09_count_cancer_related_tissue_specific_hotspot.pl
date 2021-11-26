#得../../output/cancer_total/09_count_cancer_related_tissue_specific_hotspot.txt.gz
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use Parallel::ForkManager;




 
my $fo1 = "../../output/cancer_total/09_count_cancer_related_tissue_specific_hotspot.txt.gz";
open my $O1, "| gzip >$fo1" or die $!;

print $O1 "cancer\ttissue\tnumber_of_cancer_relate_tissue_specific_hotspot\tnumber_of_cancer_relate_tissue_all_hotspot\n";


my $f1 = "../../output/cancer_total/specific/pure/tissue_specific.bed.gz";
# open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
my $f2 = "../../../GTEx/output/Tissue_total/number_of_tissue_all_hotspot.txt.gz";
open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件
my (%hash1,%hash2);

while(<$I1>)
{
    chomp;
    unless(/^h_chr/){
       my @f= split/\t/;
       my $v =join("\t",@f[0..2]);
       my $k = join("\t",@f[3,4]);
       push @{$hash1{$k}},$v;
    }
}    

while(<$I2>)
{
    chomp;
    unless(/^Tissue/){
       my @f= split/\t/;
       my $tissue = $f[0];
       my $number = $f[1];
      $hash2{$tissue}=$number;
    }
}  

foreach my $k(sort keys %hash1 ){
    my @vs = @{$hash1{$k}};
    my %hash3;
    @vs = grep { ++$hash3{$_} < 2 } @vs; 
    my $count =@vs;
    my @t =split/\t/,$k;
    my $tissue  = $t[1];
    if(exists $hash2{$tissue}){
        my $all_tissue = $hash2{$tissue};
        print $O1 "$k\t$count\t$all_tissue\n";
    }

}
