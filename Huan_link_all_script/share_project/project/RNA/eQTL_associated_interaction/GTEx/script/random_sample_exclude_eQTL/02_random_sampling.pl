#根据"/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/random_sample/01_count_number_of_eQTL_chr_in_lung_and_whole_blood.txt.gz"，对应染色体从"/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/${tissue}/${tissue}_cis_eQTL_1kg_Completion.txt.gz"随机抽取相同数目的位置作为新的QTL,并将其他位置定义为非qtl得"$output_dir/${tissue}_random_select_result_${i}.txt.gz"
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use Parallel::ForkManager;

my @tissues = ("Whole_Blood","Lung");

my $f1 = "/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/random_sample/01_count_number_of_eQTL_chr_in_lung_and_whole_blood.txt.gz";
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件

my (%hash1,%hash4);
while(<$I1>)
{
    chomp;
    unless(/^Tissue/){
        my @f = split/\t/;
        my $tissue = $f[0];
        my $chr =$f[1];
        my $number =$f[2];
        $hash1{$tissue}{$chr}=$number; #hash of hash,$hash1{$tissue}=$chr,$hash1{$tissue}为hash名, $chr为key,$number为Value
    }
}


foreach my $tissue (sort keys %hash1){
    my %hash2;
    my $f2 = "/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/${tissue}/${tissue}_cis_eQTL_1kg_Completion.txt.gz";
    # my $f2 = "/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Whole_Blood/123.txt.gz";
    open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件
    my %hash5;
    while(<$I2>)
    {
        chomp;
        unless(/^SNP_chr/){
            my @f = split/\t/;
            my $SNP_chr =$f[0];
            my $SNP_pos =$f[1];
            my $Pvalue =$f[2];
            if ($Pvalue <5e-8){#--------原本的qtl,不被random select
                push @{$hash4{$SNP_chr}},$SNP_pos;
            }
            else{
                push @{$hash2{$SNP_chr}},$SNP_pos;
            }
        }
    } 
    
    my $output_dir = "/share/data0/QTLbase/huan/GTEx/random_select_exclude_eQTL/${tissue}/random_select_result";
    if(-e $output_dir){
        print "${output_dir}\texist\n";
    }
    else{
        system "mkdir -p $output_dir";
    }
    #-----------------
    # my $pm = Parallel::ForkManager->new(10); ## 设置最大的线程数目
    for(my $i=1; $i<11;$i++){#random 取10 个看看 
        # my $pid = $pm->start and next; #开始多线程
        # my $i =1;
        print "start $tissue\t$i\n";
        my $fo1 = "$output_dir/${tissue}_random_select_result_${i}.txt.gz";
        open my $O1, "| gzip >$fo1" or die $!;
        print $O1 "SNP_chr\tSNP_pos\tPvalue\n";
        # open my $O, '>', $fo or die "$0 : failed to open output file '$fo' : $!\n";
    #---------------
        foreach my $SNP_chr(sort keys %hash2){
            if (exists $hash1{$tissue}{$SNP_chr} ){
                #------------
                my $QTL_number = $hash1{$tissue}{$SNP_chr};
                my @SNP_poss = @{$hash2{$SNP_chr}};
                #----------original qtl
                my @SNP_org_qtl = @{$hash4{$SNP_chr}};
                my %hash_u_o_q;
                @SNP_org_qtl = grep { ++$hash_u_o_q{$_} <2 } @SNP_org_qtl;
                foreach my $o_qtl(@SNP_org_qtl){
                    print $O1 "$SNP_chr\t$o_qtl\t0.05\n";
                }
                my %hash_a = map{$_=>1} @SNP_org_qtl;
                #-----------------------
                my %hash_unique;
                @SNP_poss = grep { ++$hash_unique{$_} <2 } @SNP_poss;
                my @SNP_poss_only = grep {!$hash_a{$_}} @SNP_poss; #@b_only  #-------不在@SNP_org_qtl 中的@SNP_poss
                my @new_poss = randomElem ( $QTL_number, @SNP_poss_only ) ; # pick any $QTL_number from @SNP_poss_only ，把$num和@array传递给子程序。这里是用的值传递。还有一种方式是引用传递，相当于硬链接
                # my %hash_a = map{$_=>1} @new_poss;
                # my @b_only = grep {!$hash_a{$_}} @SNP_poss; #-------不在@new_poss 中的@SNP_poss
                # my @other_snps = @b_only; #--------------
                #--------------------
                my %hash3;
                foreach my $sig_qtl(@new_poss){ #-----QTL signicifant
                    my $output = "$SNP_chr\t$sig_qtl\t5e-9";
                    print $O1 "$output\n";
                    $hash3{$sig_qtl}=1;
                }
                #--------------
                foreach my $other_snp(@SNP_poss_only){
                    unless(exists $hash3{$other_snp}){
                        my $output = "$SNP_chr\t$other_snp\t0.05";
                        print $O1 "$output\n";
                    }
                }
            }     
        }
    print "finish $tissue\t$i\n";
    # $pm->finish;  #多线程结束
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