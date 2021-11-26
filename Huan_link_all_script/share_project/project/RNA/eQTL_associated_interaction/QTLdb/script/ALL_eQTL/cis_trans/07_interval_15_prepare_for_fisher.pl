##看各种factor在#"~/project/RNA/eQTL_associated_interaction/QTLbase/output/ALL_eQTL/cis_trans/interval_15/annotation_out/01_annotation_merge_${group}_${type}.txt.gz 的fisher exact test 准备数据,
#得"/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/output/ALL_eQTL/cis_trans/interval_15/annotation_out/fisher_exact_test/prepare/07_interval_15_prepare_fisher_test.txt"
#0-based
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;

my @types= ("cis_1MB","cis_10MB","trans_1MB","trans_10MB");
# my @groups = ("hotspot","non_hotspot");

my $fo1 = "/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/output/ALL_eQTL/cis_trans/interval_15/annotation_out/fisher_exact_test/07_interval_15_prepare_fisher_test.txt";
open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
# open my $O1, "| gzip >$fo1" or die $!;
print $O1 "Type\tFactor\tNumber_of_factor_in_hotspot\tNumber_of_non_factor_in_hotspot\tNumber_of_factor_in_non_hotspot\tNumber_of_non_factor_in_non_hotspot\n";
foreach my $type(@types){
    # foreach my $group(@groups){
    my $f1 = "~/project/RNA/eQTL_associated_interaction/QTLbase/output/ALL_eQTL/cis_trans/interval_15/annotation_out/01_annotation_merge_hotspot_${type}.txt.gz";
    # open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
    open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
    #--------------------
    my $f2 = "~/project/RNA/eQTL_associated_interaction/QTLbase/output/ALL_eQTL/cis_trans/interval_15/annotation_out/01_annotation_merge_non_hotspot_${type}.txt.gz";
    # open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
    open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件
    #--------------------

    my (%hash1,%hash2,%hash3,%hash4);
    #-------hotspot
    my @factors=();
    while(<$I1>)
    {
        chomp;
        my @f = split/\t/;
        my $chr = $f[0];
        my $start = $f[1];
        my $end = $f[2];
        # my @factors = @f[3..12];
        
        if(/\bchr\b/){
        #    print $O1 "chr\tstart\tend\tscore\n";
            @factors = @f;
        }
        else{
            for (my $i=3;$i<13;$i++){
                if (abs($f[$i] -0) < 0.000000001){ # hotspot:0
                    push @{$hash1{$factors[$i]}},$f[$i];
                }
                else{# hotspot:1
                    push @{$hash2{$factors[$i]}},$f[$i];
                }
            }
        }
    }
    #----------------------------------non_hotspot
    @factors=();
    while(<$I2>)
    {
        chomp;
        my @f = split/\t/;
        my $chr = $f[0];
        my $start = $f[1];
        my $end = $f[2];
        # my @factors = @f[3..12];
        
        if(/\bchr\b/){
        #    print $O1 "chr\tstart\tend\tscore\n";
            @factors = @f;
        }
        else{
            for (my $i=3;$i<13;$i++){
                if (abs($f[$i] -0) < 0.000000001){ # non_hotspot:0
                    push @{$hash3{$factors[$i]}},$f[$i];
                }
                else{# non_hotspot:1
                    push @{$hash4{$factors[$i]}},$f[$i];
                }
            }
        }
    }

    foreach my $k(sort keys %hash1){
        if (exists $hash2{$k} && exists $hash3{$k} && exists $hash4{$k}){
            my @vs1 = @{$hash1{$k}};
            my @vs2 = @{$hash2{$k}};
            my @vs3 = @{$hash3{$k}};
            my @vs4 = @{$hash4{$k}};
            my $number1 = @vs1;
            my $number2 = @vs2;
            my $number3 = @vs3;
            my $number4= @vs4;
            print $O1 "$type\t$k\t$number1\t$number2\t$number3\t$number4\n";
        }
    }

}
