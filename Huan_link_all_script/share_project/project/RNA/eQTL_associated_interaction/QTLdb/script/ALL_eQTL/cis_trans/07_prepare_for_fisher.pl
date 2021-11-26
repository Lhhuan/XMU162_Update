##看各种factor在#"~/project/RNA/eQTL_associated_interaction/QTLdb/output/ALL_eQTL/cis_trans/fisher_exact_test/annotation_out/hotspot/05_annotation_merge_${type}_interval_${interval}_overlap_fraction_${fraction}.txt.gz" 的fisher exact test 准备数据,
#得"/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLdb/output/ALL_eQTL/cis_trans/fisher_exact_test/fisher_exact_test/interval_${interval}/07_prepare_fisher_test.txt"
#0-based
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;
use Parallel::ForkManager;

my @types= ("cis_1MB","cis_10MB","trans_1MB","trans_10MB");
# my @groups = ("hotspot","non_hotspot");
# my @intervals = (6,7,8,9,12,15,18);
my @intervals = (18,15,12,8);
# my @intervals = (18);
# my @fractions = (1,0.9,0.8,0.7,0.6,0.5);
my @fractions = (0.1);
# my @fractions = (1,0.9);
# fractions = c(0.5,0.6,0.7,0.8,0.9,1)


# my $fo1 = "/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLdb/output/ALL_eQTL/cis_trans/fisher_exact_test/annotation_out/07_prepare_fisher_test.txt";
# open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
# # open my $O1, "| gzip >$fo1" or die $!;
# print $O1 "Type\tFactor\tNumber_of_factor_in_hotspot\tNumber_of_non_factor_in_hotspot\tNumber_of_factor_in_non_hotspot\tNumber_of_non_factor_in_non_hotspot\n";
my $pm = Parallel::ForkManager->new(7); ## 设置最大的线程数目
foreach my $interval(@intervals){
    my $pid = $pm->start and next; #开始多线程
    my $fo1 = "/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLdb/output/ALL_eQTL/cis_trans/fisher_exact_test/fisher_result/interval_${interval}/07_prepare_fisher_test.txt";
    open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
    # open my $O1, "| gzip >$fo1" or die $!;
     print $O1 "Fraction\tType\tFactor\tNumber_of_factor_in_hotspot\tNumber_of_factor_in_non_hotspot\tNumber_of_non_factor_in_hotspot\tNumber_of_non_factor_in_non_hotspot\n";       
    foreach my $fraction(@fractions){
        foreach my $type(@types){
            print "$type\n";
            my $f1 = "~/project/RNA/eQTL_associated_interaction/QTLdb/output/ALL_eQTL/cis_trans/fisher_exact_test/annotation_out/hotspot/05_annotation_merge_${type}_interval_${interval}_overlap_fraction_${fraction}.txt.gz";
            # open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
            open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
            #--------------------
            my $f2 = "~/project/RNA/eQTL_associated_interaction/QTLdb/output/ALL_eQTL/cis_trans/fisher_exact_test/annotation_out/non_hotspot/05_annotation_merge_${type}_interval_${interval}_overlap_fraction_${fraction}.txt.gz";
            # open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
            open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件
            #--------------------

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
                
                # my @n1s=();
                my @n2s=();
                my @n3s=();
                my @n4s=();
                my @vs1 = @{$hash1{$k}};
                my $number1 = @vs1;
                # push @n1s,$number1;
                if (exists $hash2{$k}){
                    my @vs2 = @{$hash2{$k}};
                    my $number2 = @vs2;
                    push @n2s,$number2;
                }
                else{
                    push @n2s,0;
                }
                #-------
                if (exists $hash3{$k}){
                    my @vs3 = @{$hash3{$k}};
                    my $number3 = @vs3;
                    push @n3s,$number3;
                }
                else{
                    push @n3s,0;
                }
                #----------------
                if(exists $hash4{$k}){
                    my @vs4 = @{$hash4{$k}};
                    my $number4= @vs4;
                    push @n4s,$number4;
                }
                else{
                    push @n4s,0;
                }
                # print $O1 "$fraction\t$type\t$k\t$number1\t$n2s[0]\t$n3s[0]\t$n4s[0]\n";
                print $O1 "$fraction\t$type\t$k\t$n2s[0]\t$n4s[0]\t$number1\t$n3s[0]\n";
            }
        }
        #----------------------
    }
    $pm->finish;  #多线程结束
}