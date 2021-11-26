##看各种ROC准备数据，利用"~/project/RNA/eQTL_associated_interaction/GTEx/output/Whole_Blood/Cis_eQTL/annotation_out/interval_18/ALL/hotspot/08_annotation_merge_interval_18_overlap_cutoff_${cutoff}_chr1_22.txt.gz" "~/project/RNA/eQTL_associated_interaction/GTEx/output/Whole_Blood/Cis_eQTL/annotation_out/interval_18/ALL/non_hotspot/08_annotation_merge_interval_18_overlap_cutoff_${cutoff}_chr1_22.txt.gz"计算roc曲线的准备数据得，
#得"/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Whole_Blood/Cis_eQTL/ROC/interval_18/ALL/09_prepare_number_ROC.txt"
#0-based
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;
use Parallel::ForkManager;

my @cutoffs = ();
for (my $i=0.05;$i<0.96;$i=$i+0.05){ #对文件进行处理，把所有未定义的空格等都替换成NONE
    push @cutoffs,$i;
    # print "$i\n";
}
push @cutoffs,0.01;
push @cutoffs,0.99;

    my $fo1 = "/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Whole_Blood/Cis_eQTL/ROC/interval_18/chr1_22/09_prepare_number_ROC.txt";
    open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
    # open my $O1, "| gzip >$fo1" or die $!;
    print $O1 "cutoff\tFactor\tNumber_of_factor_in_hotspot\tNumber_of_factor_in_non_hotspot\tNumber_of_non_factor_in_hotspot\tNumber_of_non_factor_in_non_hotspot\n";    
    my $fo2 = "/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Whole_Blood/Cis_eQTL/ROC/interval_18/chr1_22/09_prepare_TPR_FPR_ROC.txt";
    open my $O2, '>', $fo2 or die "$0 : failed to open output file '$fo2' : $!\n";
    # open my $O1, "| gzip >$fo1" or die $!;
    print $O2 "cutoff\tFactor\tTPR\tFPR\n"; 
foreach my $cutoff(@cutoffs){
    my $f1 = "~/project/RNA/eQTL_associated_interaction/GTEx/output/Whole_Blood/Cis_eQTL/annotation_out/interval_18/chr1_22/hotspot/08_annotation_merge_interval_18_overlap_cutoff_${cutoff}_chr1_22.txt.gz";
    # open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
    open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
    #--------------------
    my $f2 = "~/project/RNA/eQTL_associated_interaction/GTEx/output/Whole_Blood/Cis_eQTL/annotation_out/interval_18/chr1_22/non_hotspot/08_annotation_merge_interval_18_overlap_cutoff_${cutoff}_chr1_22.txt.gz";
    # open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
    open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件
    #--------------------
    # print "$cutoff\n";
    my (%hash1,%hash2,%hash3,%hash4);
    #-------hotspot
    # print "$f1\n";
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
            for (my $i=3;$i<5;$i++){
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
            for (my $i=3;$i<5;$i++){
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
        
        # my @n1s=();
        my @n2s=();
        my @n3s=();
        my @n4s=();
        my @vs1 = @{$hash1{$k}};
        my $number1 = @vs1; #hotspot:0 FP
        # push @n1s,$number1;
        if (exists $hash2{$k}){ # hotspot:1 TP
            my @vs2 = @{$hash2{$k}};
            my $number2 = @vs2;
            push @n2s,$number2;
        }
        else{
            push @n2s,0;
        }
        #-------# non_hotspot:0 TN
        if (exists $hash3{$k}){ # non_hotspot:0
            my @vs3 = @{$hash3{$k}};
            my $number3 = @vs3;
            push @n3s,$number3;
        }
        else{
            push @n3s,0;
        }
        #----------------# non_hotspot:1 FN
        if(exists $hash4{$k}){
            my @vs4 = @{$hash4{$k}};
            my $number4= @vs4;
            push @n4s,$number4;
        }
        else{
            push @n4s,0;
        }
        my $FP = $number1;
        my $TP = $n2s[0];
        my $TN = $n3s[0];
        my $FN = $n4s[0];
        # print $O1 "$fraction\t$type\t$k\t$n2s[0]\t$n4s[0]\t$number1\t$n3s[0]\n";
        print $O1 "$cutoff\t$k\t$TP\t$FN\t$FP\t$TN\n";
        my $tpr = $TP/($TP+$FN);
        my $fpr= $FP/($FP+$TN);
        print $O2 "$cutoff\t$k\t$tpr\t$fpr\n";
    } #Number_of_factor_in_hotspot\tNumber_of_factor_in_non_hotspot\tNumber_of_non_factor_in_hotspot\tNumber_of_non_factor_in_non_hotspot
    #-----------------------------------------optimize
    # foreach my $k(sort keys %hash1){
    #     # print "$cutoff\n"; 
    #     if (exists $hash2{$k} && exists $hash3{$k} && exists $hash4{$k}){
    #         # print "$cutoff\n"; 
    #         my @vs1 = @{$hash1{$k}};
    #         my @vs2 = @{$hash2{$k}};
    #         my @vs3 = @{$hash3{$k}};
    #         my @vs4 = @{$hash4{$k}};
    #         my $FP = @vs1;# hotspot:0
    #         my $TP = @vs2;# hotspot:1
    #         my $TN = @vs3;#non_hotspot:0
    #         my $FN= @vs4;# non_hotspot:1
    #         # print "$cutoff\t$k\t$TP\t$FN\t$FP\t$TN\n";
    #         print $O1 "$cutoff\t$k\t$TP\t$FN\t$FP\t$TN\n";
    #         my $tpr = $TP/($TP+$FN);
    #         my $fpr= $FP/($FP+$TN);
    #         print $O2 "$cutoff\t$k\t$tpr\t$fpr\n";
    #     }
    #     else{
    #         print "no value\t$cutoff\n";
    #     }
    # }






#---------------------cheak发现$FN=0
#     foreach my $k(sort keys %hash1){
#         # print "$cutoff\n"; 
#         if (exists $hash2{$k} && exists $hash3{$k} ){
#             # print "$cutoff\n"; 
#             my @vs1 = @{$hash1{$k}};
#             my @vs2 = @{$hash2{$k}};
#             my @vs3 = @{$hash3{$k}};
#             # my @vs4 = @{$hash4{$k}};
#             my $FP = @vs1;# hotspot:0
#             my $TP = @vs2;# hotspot:1
#             my $TN = @vs3;#non_hotspot:0
#             my $FN= 0;# non_hotspot:1
#             # print "$cutoff\t$k\t$TP\t$FN\t$FP\t$TN\n";
#             print $O1 "$cutoff\t$k\t$TP\t$FN\t$FP\t$TN\n";
#             my $tpr = $TP/($TP+$FN);
#             my $fpr= $FP/($FP+$TN);
#             print $O2 "$cutoff\t$k\t$tpr\t$fpr\n";
#         }
#     }

}