#统计不同cutoff下hotspot的数目/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Whole_Blood/Cis_eQTL/hotspot_cis_eQTL/interval_18/，得/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Whole_Blood/Cis_eQTL/061_hotspot_in_different_cutoff_interval_18.txt
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
# my @factors = ("promoter","enhancer");
# my @groups = ("hotspot","non_hotspot");
my $fo1 = "/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Whole_Blood/Cis_eQTL/061_hotspot_in_different_cutoff_interval_18.txt";
open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
print $O1 "Cutoff\tNumber_of_hotspot\n";

foreach my $cutoff(@cutoffs){
    my $hotspot_file = "/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Whole_Blood/Cis_eQTL/hotspot_cis_eQTL/interval_18/whole_blood_segment_hotspot_cutoff_${cutoff}.bed.gz";
    my @arg1s = stat ($hotspot_file);
    my $hotspot_size = $arg1s[7]; #取$factor_hotspot_file size
    # print "$hotspot_size\n";
    # ------------------------------ factor_hotspot_size
    if ($hotspot_size >20){ #空compressed file is 20
        my $command_hotspot = "zless $hotspot_file | wc -l" ;
        my $hotspot_line_count = wc($command_hotspot);
        print $O1 "$cutoff\t$hotspot_line_count\n";
    }
    else{ #空 compressed file
        my $hotspot_line_count = 0;
        print $O1 "$cutoff\t$hotspot_line_count\n";
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

