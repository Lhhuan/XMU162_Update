#根据/home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/used_refer/${factor}.bed.gz， /home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/used_refer/non_${factor}.bed.gz，../output/Whole_Blood/Cis_eQTL/annotation/interval_18/ALL/factor/hotspot/${cutoff}/${fraction}/${factor}_whole_blood_segment_hotspot_cutoff_${cutoff}.bed.gz，../output/Whole_Blood/Cis_eQTL/annotation/interval_18/ALL/non_factor/hotspot/${cutoff}/${fraction}/non_${factor}_whole_blood_segment_hotspot_cutoff_${cutoff}.bed.gz 以factor为基数，准备计算ROC的四格表得"/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Whole_Blood/Cis_eQTL/ROC/interval_18/ALL/08_prepare_number_ROC_refine_factor_count_fraction.txt"
#extend 100bp
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
use Parallel::ForkManager;

# my @intervals = (18,15,12,9,8,7,6);
my @cutoffs = ();
for (my $i=0.05;$i<0.96;$i=$i+0.05){ #对文件进行处理，把所有未定义的空格等都替换成NONE
    push @cutoffs,$i;
    # print "$i\n";
}
push @cutoffs,0.01;
push @cutoffs,0.99;
# my @types = ("factor","non_factor");
my @factors = ("promoter","enhancer");
# my @groups = ("hotspot","non_hotspot");
my @fractions = (0.05,0.1,0.2);
my $fo1 = "/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Whole_Blood/Cis_eQTL/ROC/interval_18/ALL/08_prepare_number_ROC_refine_factor_count_fraction.txt";
open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
print $O1 "Fraction\tFactor\tCutoff\tNumber_of_factor_in_hotspot_TP\tNumber_of_factor_in_non_hotspot_FN\tNumber_of_non_factor_in_hotspot_FP\tNumber_of_non_factor_in_non_hotspot_TN\tTPR\tFPR\n";
# TP\t$FN\t$FP\t$TN
foreach my $factor(@factors){
    my $command_factor = "zless /home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/used_refer/${factor}.bed.gz | wc -l" ;
    my $command_non_factor = "zless /home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/used_refer/non_${factor}.bed.gz | wc -l" ;

    my $factor_line_count = wc($command_factor);
    my $non_factor_line_count = wc($command_non_factor);
    # print "$non_factor_line_count\t$factor\n";
    foreach my $cutoff(@cutoffs){
        foreach my $fraction(@fractions){
            my $factor_hotspot_file = "../output/Whole_Blood/Cis_eQTL/annotation/interval_18/ALL/factor/hotspot/${cutoff}/${fraction}/${factor}_whole_blood_segment_hotspot_cutoff_${cutoff}.bed.gz";
            my $non_factor_hotspot_file = "../output/Whole_Blood/Cis_eQTL/annotation/interval_18/ALL/non_factor/hotspot/${cutoff}/${fraction}/non_${factor}_whole_blood_segment_hotspot_cutoff_${cutoff}.bed.gz";
            #--------------------------------------get size of file 
            my @arg1s = stat ($factor_hotspot_file);
            my $factor_hotspot_size = $arg1s[7]; #取$factor_hotspot_file size
            my @arg2s = stat ($non_factor_hotspot_file);
            my $non_factor_hotspot_size = $arg2s[7];
            #------------------------------------
            # print "$factor_hotspot_size\t$factor_hotspot_file\n";
            # print "$non_factor_hotspot_size\t$non_factor_hotspot_file\n";
            my @factor_hotspot_line;
            my @non_factor_hotspot_line;
            #------------------------------ factor_hotspot_size
            if ($factor_hotspot_size >20){ #空compressed file is 20
                my $command_factor_hotspot = "zless $factor_hotspot_file |cut -f4,5,6| sort -u | wc -l" ;
                my $factor_hotspot_line_count = wc($command_factor_hotspot);
                push @factor_hotspot_line,$factor_hotspot_line_count;
            }
            else{ #空 compressed file
                my $factor_hotspot_line_count = 0;
                push @factor_hotspot_line,$factor_hotspot_line_count;
            }
            # #-------------------non_factor_hotspot_size
            if ($non_factor_hotspot_size >20){ #空compressed file is 20
                my $command_non_factor_hotspot = "zless $non_factor_hotspot_file |cut -f4,5,6| sort -u | wc -l" ;
                my $non_factor_hotspot_line_count = wc($command_non_factor_hotspot);
                push @non_factor_hotspot_line,$non_factor_hotspot_line_count;
            }
            else{ #空 compressed file
                my $non_factor_hotspot_line_count = 0;
                push @non_factor_hotspot_line,$non_factor_hotspot_line_count;
            }
            #--------------------------
            my $TP = $factor_hotspot_line[0];
            my $FP = $non_factor_hotspot_line[0];
            my $FN = $factor_line_count - $TP;
            my $TN =  $non_factor_line_count - $FP;
            my $tpr = $TP/($TP+$FN);
            my $fpr= $FP/($FP+$TN);
            print $O1 "$fraction\t$factor\t$cutoff\t$TP\t$FN\t$FP\t$TN\t$tpr\t$fpr\n";
        }
    }
}


sub wc{
    my $cc = $_[0]; ## 获取参数个数
    my $result = readpipe($cc);
    my @t= split/\s+/,$result;
    my $count = $t[0];
    return($count)
}
















# my @groups = ("hotspot","non_hotspot");
# my @types = ("factor","non_factor");
# my @factors = ("promoter","enhancer");

# my $pm = Parallel::ForkManager->new(10); ## 设置最大的线程数目
# foreach my $type(@types){
#     foreach my $group(@groups){
#         foreach my $cutoff(@cutoffs){

#             my $pid = $pm->start and next; #开始多线程
#             my $input_file = "../output/Whole_Blood/Cis_eQTL/${group}_cis_eQTL/interval_18/whole_blood_segment_${group}_cutoff_${cutoff}.bed.gz";
#             # my $output_file = "../../../output/ALL_${QTL}/cis_trans/interval_15/${group}/interval_${interval}_cutoff_7.3_${type}_${QTL}_segment_${group}.bed.gz";
#             my $input_file_base_name = basename($input_file);
#             # my $dir = dirname($script);
#             # print "$input_file_base_name\n";
#             my $output_dir = "../output/Whole_Blood/Cis_eQTL/annotation/interval_18/ALL/${type}/${group}/${cutoff}"; 
#             # mkdir $PMID;
#             # #------------
#             if(-e $output_dir){
#                 print "${output_dir}\texist\n";
#             }
#             else{
#                 system "mkdir -p $output_dir";
#             }
#             #------------
#             $ENV{'input_file'}  = $input_file; #设置环境变量
#             $ENV{'input_file_base_name'} = $input_file_base_name ;
#             $ENV{'output_dir'} = $output_dir ;
#             # $ENV{'fraction'} = $fraction ;
#             my $command = "bash annotation_${type}_bedtools_intersect_interval18.sh";
#             # print "$command\n";
#             system $command;
#             print "$type\t$group\t$cutoff\n";
#             $pm->finish;  #多线程结束
#         }
#     } 
# }

