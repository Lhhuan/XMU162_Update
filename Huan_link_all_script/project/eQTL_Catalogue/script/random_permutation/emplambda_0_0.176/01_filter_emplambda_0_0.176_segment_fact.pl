###过滤出"../../../../output/${tissue}/Cis_eQTL/NHP/NHPoisson_emplambda_interval_${j}_cutoff_7.3_${tissue}.txt.gz" 中emplambda$emplambda <0.176，得 "../../../../output/Tissue_merge/Cis_eQTL/hotspot_cis_eQTL/interval_18/sampling/Tissue_merge_segment_hotspot_cutoff_0_0.176.bed.gz"
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use Parallel::ForkManager;
# my @cutoffs;
my $cutoff = 0.176;
# my $cutoff = 0;
my $j = 18;
my $tissue= "Tissue_merge";

my $f1 = "/home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/NHP/NHPoisson_emplambda_interval_18_cutoff_7.3_Tissue_merge.txt.gz";
# open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件

my %hash1;
while(<$I1>)
{
    chomp;
    unless(/emplambda/){
        my @f = split/\t/;
        my $emplambda = $f[0];
        my $pos = $f[1];
        my $chr= $f[2];
        my $v= "$pos\t$emplambda";
        unless($emplambda =~/NA/){
            # if ($chr == 1){
            push @{$hash1{$chr}},$v;
            # }
        }
    }
}

print "cutoff\t$cutoff\n";
my $fo3 = "/home/huanhuan/project/eQTL_Catalogue/output/sampling/Tissue_merge_segment_hotspot_cutoff_0_0.176.bed.gz";
open my $O3, "| gzip >$fo3" or die $!;

foreach my $chr (sort keys %hash1){
    my @vs = @{$hash1{$chr}};
    my $count=0;
    my @starts=();
    my @ends=();
    my @start_line =();
    my @end_line=();
    my @results=();
    my @contents=();
    foreach my $v(@vs){
        $count++;
        # print "$count\n";
        my @t = split/\t/,$v;
        my $pos = $t[0];
        my $emplambda =$t[1];

        # my $diff = abs($emplambda - $cutoff);
        # if ($diff<0.00000000001){
        if ($emplambda <0.176){
            push @contents,$v;
            # print  "$chr\t$v\n";
            push @results,$v;
            my $result_count=@results;
            if ($result_count <2){ #第一个大于cutoff 的值

                # print "111\n";
                push @starts,$pos;
                push @ends,$pos;
                push @start_line,$count;
                push @end_line,$count;
            }                
            else{ #第二个往后大于cutoff的值
                # print "222\n";
                my $before_end_count = $end_line[0];
                my $end_diff = $count - $before_end_count;
                if ($end_diff <2 ){ #--------------和前面end行数相差1,和前面的区域相连
                    @end_line =();
                    @ends =();
                    push @end_line,$count;
                    push @ends,$pos;
                    # print "333\n";
                    # print "$count\t2222\n";
                }
                else{ #--------和前面的片段不相连,应输出前面的片段，重新定义start
                    # print $O1 "$chr\t$starts[0]\t$ends[0]\n";
                    # print "444\n";
                    @starts =();
                    @start_line =();
                    push @starts,$pos;
                    push @start_line,$count;
                    @end_line =();
                    @ends =();
                    push @end_line,$count;
                    push @ends,$pos;
                    # print "$count\t3333333\n";
                } 
            }
        }
        else{ #--------小于cutoff的位置,输出前一个hotspot
            my $content_count=@contents; #如果contents有结果，是输出的前提
            if ($content_count >0){
                my $before_end_count = $end_line[0];
                my $end_diff = $count - $before_end_count;
                if ($end_diff <2 ){ #距离hotspot 一行的位置
                    # print $O1 "$chr\t$starts[0]\t$ends[0]\n";
                    # print "555\n";
                    my $bed_ends = $ends[0]+1-1;
                    my $f_start = $starts[0]-1;
                    print $O3 "chr${chr}\t$f_start\t$bed_ends\n";
                    # print  "chr${chr}\t$starts[0]\t$bed_ends\n";
                    @contents =();#输出一次，清空@contents
                }
            }
        } 
    }
    @results=();
}
close($O3);