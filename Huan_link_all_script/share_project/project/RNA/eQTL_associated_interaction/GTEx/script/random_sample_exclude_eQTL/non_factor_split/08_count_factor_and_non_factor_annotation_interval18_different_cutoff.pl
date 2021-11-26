#根据/home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/used_refer/${factor}.bed.gz， /home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/used_refer/non_${factor}.bed.gz，"/share/data0/QTLbase/huan/GTEx/random_select_exclude_eQTL/${tissue}/annotation/interval_18/ALL/${number}/factor/${group}/${cutoff}/${factor}_${tissue}_segment_${group}_cutoff_${cutoff}.bed.gz"，"/share/data0/QTLbase/huan/GTEx/random_select_exclude_eQTL/${tissue}/annotation/interval_18/ALL/${number}/non_factor_split/${group}/${cutoff}/non_${factor}_${tissue}_segment_${group}_cutoff_${cutoff}.bed.gz" 以factor为基数，准备计算ROC的四格表得""/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/random_sample/${tissue}/Cis_eQTL/ROC/interval_18/ALL/split_non_factor/08_prepare_number_ROC_factor_count.txt"
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
push @cutoffs,0.01;
for (my $i=0.05;$i<0.7;$i=$i+0.05){ #对文件进行处理，把所有未定义的空格等都替换成NONE
    push @cutoffs,$i;
    # print "$i\n";
}
# my @types = ("factor","non_factor");
my @factors = ("promoter","enhancer","TFBS","CHROMATIN_Accessibility","HISTONE_modification","CTCF");
my @tissues = ("Whole_Blood","Lung");
my @numbers= (1..10);
my @groups = ("hotspot");
my $group = "hotspot";
my $output_dir = "/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/random_select_exclude_eQTL/Cis_eQTL/ROC/interval_18/ALL/split_non_factor/";
if(-e $output_dir){
    print "${output_dir}\texist\n";
}
else{
    system "mkdir -p $output_dir";
}
my $fo1 = "$output_dir/08_prepare_number_ROC_factor_count.txt";
open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
print $O1 "Tissue\tRandom_number\tFactor\tCutoff\tNumber_of_factor_in_hotspot_TP\tNumber_of_factor_in_non_hotspot_FN\tNumber_of_non_factor_in_hotspot_FP\tNumber_of_non_factor_in_non_hotspot_TN\tTPR\tFPR\n";
# TP\t$FN\t$FP\t$TN
foreach my $tissue(@tissues){
    foreach my $factor(@factors){
        my $command_factor = "zless /home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/used_refer/${tissue}/${factor}_union.bed.gz | wc -l" ;
        my $command_non_factor = "zless /home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/used_refer/${tissue}/non_${factor}_split_union.bed.gz | wc -l" ;

        my $factor_line_count = wc($command_factor);
        my $non_factor_line_count = wc($command_non_factor);
        # print "$non_factor_line_count\t$factor\n";
        foreach my $number(@numbers){
            foreach my $cutoff(@cutoffs){
                my $factor_hotspot_file = "/share/data0/QTLbase/huan/GTEx/random_select_exclude_eQTL/${tissue}/annotation/interval_18/ALL/${number}/factor/${group}/${cutoff}/${factor}_${tissue}_segment_${group}_cutoff_${cutoff}.bed.gz";
                my $non_factor_hotspot_file = "/share/data0/QTLbase/huan/GTEx/random_select_exclude_eQTL/${tissue}/annotation/interval_18/ALL/${number}/non_factor_split/${group}/${cutoff}/non_${factor}_${tissue}_segment_${group}_cutoff_${cutoff}.bed.gz";
                #--------------------------------------get size of file 
                my @arg1s = stat ($factor_hotspot_file);
                my $factor_hotspot_size = $arg1s[7]; #取$factor_hotspot_file size
                my @arg2s = stat ($non_factor_hotspot_file);
                my $non_factor_hotspot_size = $arg2s[7];
                #-------------------------------
                # print "$factor_hotspot_size\t$non_factor_hotspot_size\t$number\t$cutoff\t$group\n";
                # ------------------------------------
                print "$factor_hotspot_size\t$factor_hotspot_file\n";
                print "$non_factor_hotspot_size\t$non_factor_hotspot_file\n";
                my @factor_hotspot_line;
                my @non_factor_hotspot_line;
                #-------------------------------factor_hotspot_size
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
                print $O1 "$tissue\t$number\t$factor\t$cutoff\t$TP\t$FN\t$FP\t$TN\t$tpr\t$fpr\n";

            }
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
