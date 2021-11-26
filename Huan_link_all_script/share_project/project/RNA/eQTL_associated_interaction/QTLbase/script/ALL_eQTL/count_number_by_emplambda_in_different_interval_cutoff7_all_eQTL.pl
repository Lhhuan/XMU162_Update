##统计../output/ALL_eQTL/NHPoisson_emplambda_interval_${i}_all_eQTLbase.txt中 中每个阶段eqtl的数目，
#得../output/ALL_eQTL/count_number_by_emplambda_in_different_interval_cutoff_7_all_eQTL.txt
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/sum/;
use List::Util qw/max min/;

my %hash1;


my $fo1 = "../../output/ALL_eQTL/count_number_by_emplambda_in_different_interval_cutoff_7_all_eQTL.txt";
open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
my $header = "emplambda\tnumber\tinterval";
print $O1 "$header\n";

my @intervals;
for (my $i=2000;$i<9000;$i=$i+2000){
    push @intervals, $i;
}
for (my $i=10000;$i<100001;$i=$i+10000){
    push @intervals, $i;
}

foreach my $i(@intervals) { #interval
    my $f1 = "../../output/ALL_eQTL/NHPoisson_emplambda_interval_${i}_all_eQTLbase.txt";
    open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n";
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
        open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n";
        # print "$j\n";
        while(<$I1>)
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
