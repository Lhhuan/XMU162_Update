###"/share/data0/QTLbase/huan/GTEx/random_select/${tissue}/random_select_result/${tissue}_random_select_result_${number}.txt.gz" 时得不同cutoff下的hotspot(segment),得"/share/data0/QTLbase/huan/GTEx/random_select/${tissue}/hotspot/${number}/interval_18/${tissue}_segment_hotspot_cutoff_${cutoff}.bed.gz"
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use Parallel::ForkManager;
my @cutoffs;
# push @cutoffs,0.01;
for (my $i=0.1;$i<0.31;$i=$i+0.01){ #对文件进行处理，把所有未定义的空格等都替换成NONE
    push @cutoffs,$i;
    # print "$i\n";
}

# push @cutoffs,0.99;
# my @cutoffs = (0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7)

my @numbers = (1..10);
my @tissues = ("Whole_Blood","Lung");
# my @interval = (18);
foreach my $tissue(@tissues){
    foreach my $number(@numbers){
        print "$tissue\t$number\n";
        # my $f1 = "../output/Whole_Blood/Cis_eQTL/NHP/interval_${j}_chr1.txt.gz";
        # NHPoisson_emplambda_interval_18_cutoff_7.3_Whole_Blood_1.txt
        my $f1 = "/share/data0/QTLbase/huan/GTEx/random_select/${tissue}/NHP/NHPoisson_emplambda_interval_18_cutoff_7.3_${tissue}_${number}.txt.gz";
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

        my $output_dir = "/share/data0/QTLbase/huan/GTEx/random_select/${tissue}/hotspot/${number}/interval_18";
        # #------------
        if(-e $output_dir){
            print "${output_dir}\texist\n";
        }
        else{
            system "mkdir -p $output_dir";
        }
        #--------------------
        my $pm = Parallel::ForkManager->new(20);
        foreach my $cutoff(@cutoffs){
            print "cutoff\t$cutoff\n";
            my $pid = $pm->start and next; #开始多线程
            my $fo3 = "$output_dir/${tissue}_segment_hotspot_cutoff_${cutoff}.bed.gz";
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
                    if ($emplambda >= $cutoff){
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
                                my $bed_ends = $ends[0]+1;
                                print $O3 "chr${chr}\t$starts[0]\t$bed_ends\n";
                                # print  "chr${chr}\t$starts[0]\t$bed_ends\n";
                                @contents =();#输出一次，清空@contents
                            }
                        }
                    } 
                }
                @results=();
            }
            close($O3);
            $pm->finish;  #多线程结束
        }
    }
}