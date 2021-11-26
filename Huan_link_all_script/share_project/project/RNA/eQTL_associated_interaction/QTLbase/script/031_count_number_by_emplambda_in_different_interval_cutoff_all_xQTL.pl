###统计各种QTL../output/ALL_${xQTL}/QTLbase_NHPoisson_emplambda_interval_${i}_cutoff_${cutoff}_all_${xQTL}.txt.gz中 中每个阶段eqtl的数目，
#得../output/ALL_${xQTL}/count_number_by_emplambda_in_different_interval_${cutoff}_all_${xQTL}.txt
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/sum/;
use List::Util qw/max min/;


my @QTLs = ("eQTL","caQTL","edQTL","hQTL","mQTL","pQTL","reQTL","sQTL","QTL","cerQTL","lncRNAQTL","metaQTL","miQTL","riboQTL");
# my @cutoffs=("7","24.25");
my @cutoffs=("7","7.3");

my @intervals;
for (my $i=1000;$i<5001;$i=$i+1000){
    push @intervals, $i;
}
foreach my $xQTL(@QTLs){
    print "$xQTL\n";
    foreach my $cutoff(@cutoffs){
        
        my $fo1 = "../output/ALL_${xQTL}/count_number_by_emplambda_in_different_interval_${cutoff}_all_${xQTL}.txt";
        open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
        my $header = "emplambda\tnumber\tinterval";
        print $O1 "$header\n";
        foreach my $i(@intervals) { #interval
            my %hash1;
            my $f1 = "../output/ALL_${xQTL}/QTLbase_NHPoisson_emplambda_interval_${i}_cutoff_${cutoff}_all_${xQTL}.txt.gz";
            open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
            # open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n";
            my @values=();
            while(<$I1>)
            {
                chomp;
                unless(/^emplambda/){
                    my @f=split/\t/;
                    my $emplambda = $f[0];
                    my $t =$f[1];
                    my $chr =$f[2];
                    unless($emplambda =~/NA/){
                        push @values,$emplambda;
                    }
                }
            }
            my $max_emplambda =max @values; #find the max emplambda
            print "$max_emplambda\n";

            #_------------------------------
            my $score_interval = 0.01;
            my $range = $max_emplambda+$score_interval;
            for (my $j=$score_interval;$j<$range;$j=$j+$score_interval){ #emplambda 间隔
                my $f2 =$f1;
                open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件
                # print "$j\n";
                while(<$I2>)
                {
                    chomp;
                    unless(/^emplambda/){
                        my @f=split/\t/;
                        my $emplambda = $f[0];
                        my $t =$f[1];
                        my $chr =$f[2];
                        unless($emplambda =~/NA/){
                            my $start = $j-$score_interval;
                            my $end =$j;
                            my $v = "$chr\t$t";
                            if ($emplambda>$start && $emplambda <= $end){
                                push @{$hash1{$end}},$v;
                                
                            }
                        }
                    }
                }
                
            }

            #----------------------
            foreach my $k(sort keys %hash1){
                my @vs = @{$hash1{$k}};
                my $number= @vs;
                print $O1 "$k\t$number\t$i\n";
            }
        }

    }
}