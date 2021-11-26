#将../output/tissue_pop_used_to_clump/$used_to_clump 和../output/tissue_pop_used_to_clump/${used_to_clump}.clump merge 起来，并且计算每个LD的长度得../output/clump_region/$used_to_clump
#得总文件../output/all_qtl_clump_locus.txt.gz

#!/usr/bin/perl
use warnings;
use strict; 
use utf8;

my $r_square = "0.3";
my $fo1 = "../output/all_qtl_clump_locus_r_square${r_square}.txt.gz";
open my $O1, "| gzip >$fo1" or die $!;
print $O1 "SNP\tMapped_gene\tTrait_chr\tTrait_start\tTrait_end\tP\tSourceid\tTag_snp\tLocus_region\tLocus_distance\tQTL_type\tTissue\tPopulation\n";

my $fo3 = "../output/unique_by_qtl_pop_locus_r_square${r_square}.txt.gz";
open my $O3, "| gzip >$fo3" or die $!;
print $O3 "Locus_region\tLocus_distance\tQTL_type\tPopulation\n";
# open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";   
#-------------------------------------------------------------#获取
my $dir = "../output/tissue_pop_used_to_clump/";
opendir (DIR, $dir) or die "can't open the directory!";
my @files = readdir DIR; #获取一个文件夹下的所有文件
my @tissues;
# my $suffix = "QTL.txt.gz"; #文件名后缀QTL.txt
foreach my $file(@files){ #分QTL进行
    if($file =~/QTL/){
        my $f1 = "../output/tissue_pop_clump_${r_square}/${file}.clumped";
        if (-e $f1){ #因为有些文件没有对应的clumped文件存在，所以首先判断$f1存在，
            my $f2 = "../output/tissue_pop_used_to_clump/$file";
            open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n";
            open my $I2, '<', $f2 or die "$0 : failed to open input file '$f2' : $!\n";  

            my $fo2 = "../output/clump_region_${r_square}/${file}.txt.gz";  
            open my $O2, "| gzip >$fo2" or die $!;

            print $O2 "SNP\tMapped_gene\tTrait_chr\tTrait_start\tTrait_end\tP\tSourceid\tTag_snp\tLocus_region\tLocus_distance\n";
            my %hash1;
            my %hash3;

            while(<$I1>)
            {
                chomp;
                unless(/SNP/){
                    my @f =split/\s+/;
                    for (my $i=0;$i<13;$i++){ #对文件进行处理，把所有未定义的空格等都替换成NONE
                        unless(defined $f[$i]){
                        $f[$i] = "NA";
                        }
                    unless($f[$i]=~/\w/){$f[$i]="NA"} #对文件进行处理，把所有定义的没有字符的都替换成NULL
                    }
                    my $tag_snp = $f[3];
                    unless ($tag_snp =~/NA/){ #跳过空行
                        my @tags = split/\_/,$tag_snp;
                        my $chr = $tags[0];#ld 发生在同一条染色体上，
                        my $snp2 = $f[-1];
                        my @ts= split/\,/,$snp2;
                        
                        push @ts,$tag_snp; #为了方便计算，把tag snp也 push进数组
                        my %hash2;
                        @ts = grep { ++$hash2{$_} < 2 } @ts; #数组内去重

                        my @all_unique_snp;
                        foreach my $t(@ts){ #此处为了对（）进行替换，因为不确定是(1),所以用下面策略
                            unless($t =~/NONE/){
                                $t =~s/\(.*+//g;
                                my @p = split/\_/,$t;
                                my $pos = $p[1];
                                push @all_unique_snp,$pos;
                            }
                        }
                        #-----------------------------
                        my $number = @all_unique_snp;
                        if ($number >1){ #该经过 clump,tag snp对应的locus里面有被clump的snp
                            my @sorted_all_unique_snp = sort { $a <=> $b } @all_unique_snp; #对数组内的数字进行降序排列
                            my $low = $sorted_all_unique_snp[0];
                            my $up = $sorted_all_unique_snp[-1];
                            my $distance = $up - $low;

                            my $locus_region = "${chr}_${low}_${up}";
                            foreach my $snp (@sorted_all_unique_snp){
                                my $snp_in_locus = "${chr}_${snp}"; #还原snp的格式
                                my $v= "$tag_snp\t$locus_region\t$distance";
                                $hash1{$snp_in_locus}=$v;
                            }
                        }
                        else{ #locus里面只有tag snp,@sorted_all_unique_snp 只有一个元素
                            my $distance =0;
                            my $low = $all_unique_snp[0];
                            my $up = $all_unique_snp[0];
                            my $locus_region = "${chr}_${low}_${up}";
                            my $snp = $all_unique_snp[0];
                            my $snp_in_locus = "${chr}_${snp}"; #还原snp的格式
                            my $v= "$tag_snp\t$locus_region\t$distance";
                            my $k= $snp_in_locus;
                            $hash1{$snp_in_locus}=$v;
                        }
                    }
                }
            }
#——------------------------------------
            my @info= split/\_/,$file; # exact QTL_type, tissue, population
            my $info_number = @info;
            if ($info_number >3){ #tissue中间有_ 比如"Large_Intestine-Colon"
                my $qtl_type = $info[0];
                my $number_for_tissue = $info_number-2; #获取数组倒数第二个元素
                my $tissue =join("_",@info[1..$number_for_tissue]);
                my $population = $info[-1];
                while(<$I2>)
                {
                    chomp;
                    unless(/^SNP/){
                        my @f =split/\t/;
                        my $snp = $f[0];
                        if (exists $hash1{$snp}){
                            my $v = $hash1{$snp};
                            my @t= split/\t/,$v;
                            my $tag_snp =$t[0];
                            my $locus_region = $t[1];
                            my $distance =$t[2];
                            my $output3 = "$locus_region\t$distance\t$qtl_type\t$population";
                            print $O1 "$_\t$v\t$qtl_type\t$tissue\t$population\n";
                            print $O2 "$_\t$v\n";
                            unless (exists $hash3{$output3}){
                                $hash3{$output3}=1;
                                print $O3 "$output3\n";
                            }
                        }
                    }
                }
            }
            else{ #文件名是三个元素组成
                my $qtl_type = $info[0];
                my $tissue = $info[1];
                my $population = $info[-1];
                while(<$I2>)
                {
                    chomp;
                    unless(/^SNP/){
                        my @f =split/\t/;
                        my $snp = $f[0];
                        if (exists $hash1{$snp}){
                            my $v = $hash1{$snp};
                            my @t= split/\t/,$v;
                            my $tag_snp =$t[0];
                            my $locus_region = $t[1];
                            my $distance =$t[2];
                            my $output3 = "$locus_region\t$distance\t$qtl_type\t$population";
                            print $O1 "$_\t$v\t$qtl_type\t$tissue\t$population\n";
                            print $O2 "$_\t$v\n";
                            unless (exists $hash3{$output3}){
                                $hash3{$output3}=1;
                                print $O3 "$output3\n";
                            }
                        }
                    }
                }
            }
            print "finish $f1\n";
        }
    }    
}