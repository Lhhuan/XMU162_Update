#根据../../output/random_segmention/01_count_number_and_length_of_hotspot_chr_in_lung_and_whole_blood.txt.gz 从"/home/huanhuan/ref_data/UCSC/hg19.chrom1_22.sizes"中随机抽取染色体特异性特定长度的片段
#得
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use Parallel::ForkManager;

# my @tissues = ("Whole_Blood","Lung");
my @tissues = ("Lung");
my %hash1;
my $f1 = "/home/huanhuan/ref_data/UCSC/hg19.chrom1_22.sizes";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 

while(<$I1>)
{
    chomp;
    unless(/^Tissue/){
        my @f = split/\t/;
        my $chr = $f[0];
        # $chr =~ s/chr//g;
        my $size =$f[1];
        $hash1{$chr}=$size;
    }
}

my @cutoffs;
push @cutoffs,0.01;
# for (my $i=0.05;$i<0.7;$i=$i+0.05){ #对文件进行处理，把所有未定义的空格等都替换成NONE
#     push @cutoffs,$i;
#     # print "$i\n";
# }

foreach my $tissue(@tissues){
    foreach my $cutoff(@cutoffs){
        my %hash2;
        my $f2 = "../../output/random_segmention/01_count_number_and_length_of_hotspot_chr_in_${tissue}.txt.gz";
        open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件

        while(<$I2>)
        {
            chomp;
            unless(/^Chr/){
                my @f = split/\t/;
                unless(/^Chr/){
                    my $chr =$f[0];
                    my $Length =$f[1];
                    my $Number =$f[2];
                    my $Cutoff =$f[3];
                    if (abs($Cutoff-$cutoff)<0.00000001){
                        $hash2{$chr}{$Length}=$Number;
                        # print "$Cutoff\t$cutoff\n";
                    }
                    
                }
            }
        }

        foreach my $chr(sort keys %hash1){
            my $size = $hash1{$chr};
            my @chr_array =(1..$size);
            if (exists $hash2{$chr}){
                my @Lengths=();
                my @Numbers=();
                foreach my $length( keys %{$hash2{$chr}}){
                    my $num = $hash2{$chr}{$length};
                    push @Lengths,$length;
                    push @Numbers,$num;
                }
                my $number = sum @Numbers;
                my $all_number = $number +1; #总长度比要选的点多1，防止最后一个点不能作为start
                # print "$number\t$all_number\n";
                my @random_starts = randomElem ( $all_number, @chr_array ) ; # pick any $QTL_number from @SNP_poss ，把$num和@array传递给子程序。这里是用的值传递。还有一种方式是引用传递，相当于硬链接
                my @sorted_random_starts = sort {$a <=> $b} @random_starts; #对数组内元素按照数字升序排序
                for (my $i=1,$i<$all_number,$i++){
                    my @
                    if ()
                }
            }

        }





    }
}



sub randomElem { #随机取
    my ($want, @array) = @_ ;
    my (%seen, @ret);
    while ( @ret != $want ) {
        my $num = abs(int(rand(@array))); #@array 是指数组的长度，而$#array是指最后一个索引，由于rand的特殊性，如果用$#array会导致取不到最后一个值。
        if ( ! $seen{$num} ) { 
            ++$seen{$num};
            push @ret, $array[$num];
        }
    }
    return @ret;     
}







# foreach my $tissue (sort keys %hash1){
#     my %hash2;
#     my $f2 = "/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/${tissue}/${tissue}_cis_eQTL_1kg_Completion.txt.gz";
#     # my $f2 = "/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Whole_Blood/123.txt.gz";
#     open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件

#     while(<$I2>)
#     {
#         chomp;
#         unless(/^SNP_chr/){
#             my @f = split/\t/;
#             my $SNP_chr =$f[0];
#             my $SNP_pos =$f[1];
#             my $Pvalue =$f[2];
#             push @{$hash2{$SNP_chr}},$SNP_pos;

            
#         }
#     } 
    
#     my $output_dir = "/share/data0/QTLbase/huan/GTEx/random_select/${tissue}/random_select_result";
#     if(-e $output_dir){
#         print "${output_dir}\texist\n";
#     }
#     else{
#         system "mkdir -p $output_dir";
#     }
#     #-----------------
#     # my $pm = Parallel::ForkManager->new(10); ## 设置最大的线程数目
#     for(my $i=1; $i<11;$i++){#random 取10 个看看 
#         # my $pid = $pm->start and next; #开始多线程
#         # my $i =1;
#         print "start $tissue\t$i\n";
#         my $fo1 = "$output_dir/${tissue}_random_select_result_${i}.txt.gz";
#         open my $O1, "| gzip >$fo1" or die $!;
#         print $O1 "SNP_chr\tSNP_pos\tPvalue\n";
#         # open my $O, '>', $fo or die "$0 : failed to open output file '$fo' : $!\n";
#     #---------------
#         foreach my $SNP_chr(sort keys %hash2){
#             if (exists $hash1{$tissue}{$SNP_chr} ){
#                 #------------
#                 my $QTL_number = $hash1{$tissue}{$SNP_chr};
#                 my @SNP_poss = @{$hash2{$SNP_chr}};
#                 my %hash_unique;
#                 @SNP_poss = grep { ++$hash_unique{$_} <2 } @SNP_poss;
#                 my @new_poss = randomElem ( $QTL_number, @SNP_poss ) ; # pick any $QTL_number from @SNP_poss ，把$num和@array传递给子程序。这里是用的值传递。还有一种方式是引用传递，相当于硬链接
#                 # my %hash_a = map{$_=>1} @new_poss;
#                 # my @b_only = grep {!$hash_a{$_}} @SNP_poss; #-------不在@new_poss 中的@SNP_poss
#                 # my @other_snps = @b_only; #--------------
#                 #--------------------
#                 my %hash3;
#                 foreach my $sig_qtl(@new_poss){ #-----QTL signicifant
#                     my $output = "$SNP_chr\t$sig_qtl\t5e-9";
#                     print $O1 "$output\n";
#                     $hash3{$sig_qtl}=1;
#                 }
#                 #--------------
#                 foreach my $other_snp(@SNP_poss){
#                     unless(exists $hash3{$other_snp}){
#                         my $output = "$SNP_chr\t$other_snp\t0.05";
#                         print $O1 "$output\n";
#                     }
#                 }
#             }     
#         }
#     print "finish $tissue\t$i\n";
#     # $pm->finish;  #多线程结束
#     }
# }


# sub randomElem { #随机取
#     my ($want, @array) = @_ ;
#     my (%seen, @ret);
#      while ( @ret != $want ) {
#      my $num = abs(int(rand(@array))); #@array 是指数组的长度，而$#array是指最后一个索引，由于rand的特殊性，如果用$#array会导致取不到最后一个值。
#       if ( ! $seen{$num} ) { 
#         ++$seen{$num};
#         push @ret, $array[$num];
#       }
#      }
#     return @ret;     
# }