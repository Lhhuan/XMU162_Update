#统计lung和 whole_blood ${dir}/${tissue}${suffix}中每条染色体的eQTL的数目，得../../output/random_sample/01_count_number_of_eQTL_chr_in_lung_and_whole_blood.txt.gz
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;

my $fo1 = "../../output/random_sample/01_count_number_of_eQTL_chr_in_lung_and_whole_blood.txt.gz";
open my $O1, "| gzip >$fo1" or die $!;

print $O1 "Tissue\tChr\tnumber\n";

my %hash1;

for (my $i=1;$i<23;$i++){
    my $k = $i;
    $hash1{$k}=1;
}

my $dir = "../../data/GTEx_Analysis_v8_eQTL_hg19";
my $suffix = ".v8.signif_variant_gene_pairs.txt.gz";
my @tissues = ("Whole_Blood", "Lung");

foreach my $tissue (@tissues){
    my $f1 ="${dir}/${tissue}${suffix}";
    open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件

    my(%hash2,%hash3,%hash4);
    while(<$I1>){
        chomp;
        my @f = split/\t/;
        unless(/^variant_id/){
            my $variant_id = $f[0];
            my $chr = $f[1];
            my $pos = $f[2];
            my $Pvalue =$f[-6];
            if($Pvalue <5e-8){ #----significant qtl
                if(exists $hash1{$chr}){
                    push @{$hash2{$chr}},$pos;
                }
            }

        }
    }

    foreach my $k(sort keys %hash2){
        my @poss = @{$hash2{$k}};
        my %hash; 
        @poss = grep { ++$hash{$_}< 2 } @poss;
        my $number = @poss;
        print $O1 "$tissue\t$k\t$number\n";
    }
}