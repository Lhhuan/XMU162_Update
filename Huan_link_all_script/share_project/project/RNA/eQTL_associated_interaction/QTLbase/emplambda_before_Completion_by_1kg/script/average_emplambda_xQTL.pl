#将../output/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_eQTL_caQTL_edQTL_hQTL_mQTL_pQTL_reQTL_sQTL_allQTL.txt.gz 按照特定距离将emplambda取平均
#得../output/average_emplambda_xQTL.txt.gz
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;

my $f1 = "../output/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_eQTL_caQTL_edQTL_hQTL_mQTL_pQTL_reQTL_sQTL_allQTL.txt.gz";
# open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件

#------------------------------------

my $fo1 = "../output/emplambda_in_chr_min_max_position.txt";
open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
my $fo2 = "../output/emplambda_in_chr_max_position.txt";
open my $O2, '>', $fo2 or die "$0 : failed to open output file '$fo2' : $!\n";

print $O1 "chr\tmin\tmax\n";

my (%hash1,%hash2);
while(<$I1>)
{
    chomp;
    unless(/emplambda/){
        my @f = split/\t/;
        my $emplambda =$f[0];
        my $t =$f[1]; 
        my $pos = $t; #snp pos
        $pos =~ s/"//g;
        my $chr =$f[2]; #snp chr
        my $xQTL = $f[3];
        my $k1 = $chr;
        unless ($emplambda =~/NA/){ 
            push @{$hash1{$chr}},$pos;
        }
    }
}

my %hash3;
foreach my $k1(sort keys %hash1){
    my @pos = @{$hash1{$k1}};
    my $max_pos = max @pos;
    my $min_pos =min @pos;
    my $v = "$min_pos\t$max_pos";
    print $O1 "$k1\t$v\n";
    print $O2 "$k1\t$$max_pos\n";
}