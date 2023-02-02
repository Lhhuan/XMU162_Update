#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
use List::MoreUtils ':all';
# use Parallel::ForkManager;


my @states=("1_TssA","2_TssAFlnk","3_TxFlnk","4_Tx","5_TxWk","6_EnhG","7_Enh","8_ZNF/Rpts","9_Het","10_TssBiv","11_BivFlnk","12_EnhBiv","13_ReprPC","14_ReprPCWk","15_Quies");

foreach my $state(@states){
    my $state1 = $state;
    $state1 =~s/\//_/g;
    my $f1 = "./all_hg38/E062_15_coreMarks_hg38lift_dense.bed.gz";
    open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); 
    my $fo2 = "./EO62/E062_15_coreMarks_hg38lift_dense_${state1}.bed.gz";
    open my $O2, "| gzip >$fo2" or die $!;
    my %hash1;

    while(<$I1>)
    {
        chomp;
        unless(/^track/){
            my @f= split/\t/;
            my $st = $f[3];
            if($st eq $state){
                # my $output=join("\t",@f[0..3]);
                my $output=join("\t",@f[0..2]);
                print $O2 "$output\n";
            }
        }
    }
    close($O2);
}


