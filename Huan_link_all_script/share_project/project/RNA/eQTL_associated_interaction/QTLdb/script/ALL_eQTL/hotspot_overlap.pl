#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
# use Parallel::ForkManager;

my @interval = (6,7,8,9,12,15,18,19);
my (%hash6,%hash7,%hash8,%hash9,%hash12,%hash15,%hash19);
my $fo1 = "../../output/ALL_eQTL/hotspot.txt.gz";
open my $O1, "| gzip >$fo1" or die $!;
print $O1 "chr\tpos\n";

foreach my $j(@interval){
    # print "%hash${j}\n";
    my $f1 = "../../output/ALL_eQTL/NHPoisson_emplambda_interval_${j}cutoff_7.3_all_eQTL.txt.gz";
    # open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
    open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
    my $xQTL = "eQTL" ;
    while(<$I1>)
    {
        chomp;
        unless(/emplambda/){
            my @f = split/\t/;
            my $emplambda = $f[0];
            my $pos = $f[1];
            my $chr= $f[2];
            my $k = "$chr\t$pos";
            unless($emplambda =~/NA/){
                # if ($j == 6 || $j ==7 || $j=8 ||$j==9 ||$j=12 || $j=15){
                if ($j == 6 && $emplambda >=0.2){
                    $hash6{$k} =1;
                }
                elsif($j == 7 && $emplambda >=0.2){
                    $hash7{$k} =1;
                }
                elsif($j == 8 && $emplambda >=0.2){
                    $hash8{$k} =1;
                }
                elsif($j == 9 && $emplambda >=0.2){
                    $hash9{$k} =1;
                }
                elsif($j == 12 && $emplambda >=0.2){
                    $hash12{$k} =1;
                }
                elsif($j == 15 && $emplambda >=0.2){
                    $hash15{$k} =1;
                }
                else{ #j=19
                    if ($emplambda >=0.4){
                        $hash19{$k} =1;
                    }
                }
            }
        }
    }
}

foreach my $k(sort keys %hash6){
    if (exists $hash7{$k} && $hash8{$k} &&$hash9{$k} &&$hash12{$k} && $hash15{$k} &&$hash19{$k}){
        print $O1 "$k\n";
    }
}
