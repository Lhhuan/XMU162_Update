#interval_18 时，对."/share/data0/QTLbase/huan/GTEx/random_select/${tissue}/${group}/${number}/interval_18/${tissue}_segment_${group}_cutoff_${cutoff}.bed.gz"用annotation_${type}_${tissue}_bedtools_intersect_interval18.sh进行annotation,得${tissue}_segment_${group}_cutoff_${cutoff}.bed.gz
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
# push @cutoffs,0.01;
for (my $i=0.1;$i<0.31;$i=$i+0.01){ #对文件进行处理，把所有未定义的空格等都替换成NONE
    push @cutoffs,$i;
    # print "$i\n";
}
# push @cutoffs,0.99;
my @groups = ("hotspot");
my @types = ("non_factor_split","factor");
my @tissues = ("Whole_Blood","Lung");
my @numbers= (1..10);


my $pm = Parallel::ForkManager->new(10); ## 设置最大的线程数目
foreach my $tissue(@tissues){
    foreach my $type(@types){
        foreach my $group(@groups){
            foreach my $number(@numbers){
                foreach my $cutoff(@cutoffs){
                    my $pid = $pm->start and next; #开始多线程
                    my $input_file = "/share/data0/QTLbase/huan/GTEx/random_select/${tissue}/${group}/${number}/interval_18_filter/6_5000/${tissue}_segment_${group}_cutoff_${cutoff}.bed.gz";
                    # my $output_file = "../../../output/ALL_${QTL}/cis_trans/interval_15/${group}/interval_${interval}_cutoff_7.3_${type}_${QTL}_segment_${group}.bed.gz";
                    my $input_file_base_name = basename($input_file);
                    # my $dir = dirname($script);
                    # print "$input_file_base_name\n";
                    my $output_dir = "/share/data0/QTLbase/huan/GTEx/random_select/${tissue}/annotation/interval_18_filter/6_5000/ALL/${number}/${type}/${group}/${cutoff}"; 
                    # mkdir $PMID;
                    # #------------
                    if(-e $output_dir){
                        print "${output_dir}\texist\n";
                    }
                    else{
                        system "mkdir -p $output_dir";
                    }
                    #------------
                    $ENV{'input_file'}  = $input_file; #设置环境变量
                    $ENV{'input_file_base_name'} = $input_file_base_name ;
                    $ENV{'output_dir'} = $output_dir ;
                    # $ENV{'fraction'} = $fraction ;
                    my $command = "bash annotation_${type}_${tissue}_bedtools_intersect_interval18.sh";
                    # print "$command\n";
                    system $command;
                    print "$type\t$group\t$cutoff\n";
                    $pm->finish;  #多线程结束
                }
            } 
        }
    }
}

