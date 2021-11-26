#根据"../../output/random_segmention/03_count_average_length_and_all_number_of_hotspot_chr_in_${tissue}.txt.gz"和"/home/huanhuan/ref_data/UCSC/hg19.chrom1_22.sizes" split 染色体并random取hotspot的均长，得../


#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/sum/;
use POSIX;
use List::Util qw( min max );

my $f1 = "/home/huanhuan/ref_data/UCSC/hg19.chrom1_22.sizes";
# open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 

my %hash1;
while(<$I1>){
    chomp;
    my @f = split/\t/;
    my $chr = $f[0];
    my $size = $f[1];
    $hash1{$chr}=$size;
}


# my @tissues = ("Lung");
my @tissues = ("Lung","Whole_Blood");
my @cutoffs;
# push @cutoffs,0.1;

# for (my $i=0.05;$i<0.21;$i=$i+0.05){ #对文件进行处理，把所有未定义的空格等都替换成NONE
#     push @cutoffs,$i;
#     # print "$i\n";
# }


for (my $i=0.1;$i<0.7;$i=$i+0.05){ #对文件进行处理，把所有未定义的空格等都替换成NONE
    push @cutoffs,$i;
    # print "$i\n";
}



foreach my $tissue(@tissues){
    my (%hash2,%hash3);
    
    my $f2 = "../../../output/random_segmention/03_count_average_length_and_all_number_of_hotspot_chr_in_${tissue}.txt.gz";
    open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件
    while(<$I2>){
        chomp;
        unless(/^Cutoff/){
            my @f = split/\t/;
            my $Cutoff =$f[0];
            my $Chr = $f[1];
            my $Average_Length = $f[2];
            my $hotspot_Number = $f[-1];
            # if( $Chr =~/\bchr1\b/ ||$Chr =~/\bchr22\b/ ){
            # if( $Chr =~/\bchr10\b/){
                my $k = "$Chr\t$Cutoff";
                my $v = "$Average_Length\t$hotspot_Number";
                # print "$k\t$v\n";
                $hash2{$k}=$v;
                # print "$hotspot_Number\n";
                if (exists $hash1{$Chr}){
                    my $chr_size = $hash1{$Chr};
                    # my $fo1 = "./tmp/chr_size.txt";
                    # # open my $O1, "| gzip >$fo1" or die $!;
                    # open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
                    # print $O1 "$Chr\t$chr_size\n";
                    # my $command2 = "bedtools makewindows -g  ./tmp/chr_size.txt -n $hotspot_Number | gzip > /share/data0/QTLbase/huan/GTEx/random_segment/mean_hotspot/split_chr_size/split_${Cutoff}_${Chr}_${tissue}.bed.gz";
                    # system $command2;
                }
            # }
        }
    }
    foreach my $Cutoff(@cutoffs){
        # print "$Cutoff\n";
    #----------------------
        # for(my $i=1; $i<11;$i++){
            my $i=6;
            my $out_dir = "/share/data0/QTLbase/huan/GTEx/random_segment/mean_hotspot/${tissue}/hotspot/${i}/interval_18/";
            if(-e $out_dir){
                # print "${out_dir}\texist\n";
            }
            else{
                system "mkdir -p $out_dir";
            }
            my $out_file = "${tissue}_segment_hotspot_cutoff_${Cutoff}.bed.gz";
            my $fo2 = "${out_dir}/${out_file}";
            # open my $O1, "| gzip >$fo1" or die $!;
            open my $O2, '>', $fo2 or die "$0 : failed to open output file '$fo2' : $!\n";
            foreach my $k(sort keys %hash2){
                my @t= split/\t/,$k;
                my $Chr =$t[0];
                my $cutoff =$t[1];
                my $v= $hash2{$k};
                my @ss = split/\t/,$v;
                my $Average_Length = $ss[0];
                my $hotspot_Number = $ss[1];
                # my $kk3 = "$Cutoff\t$Chr\t$tissue";
                if (abs($cutoff -$Cutoff)<0.0000001){  #不同chr,有些cutoff是不存在的
                # if ($cutoff == $Cutoff){
                print "$k\n";
                    print "start\t$i\t$Chr\t$tissue\t$Cutoff\n";
                    my $f2 = "/share/data0/QTLbase/huan/GTEx/random_segment/mean_hotspot/split_chr_size/split_${Cutoff}_${Chr}_${tissue}.bed.gz";
                    open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件
                    # -----------
                    while(<$I2>){
                        chomp;
                        my @f = split/\t/;
                        my $chr_seg = $f[0];
                        my $start = $f[1];
                        my $end = $f[2];
                        my $allow_end_to_random = $end - $Average_Length-1;
                        my @ramdom_nums  =($start..$allow_end_to_random);
                        my $random_count = 2;
                        my @random_starts = randomElem ($random_count, @ramdom_nums);
                        # my $random_start = randomElem ($random_count, @ramdom_nums);
                        # print  "$random_starts[0]\n";
                        # print "$random_start\n";
                        my $random_start = $random_starts[0];
                        my $random_end = $random_start +$Average_Length;
                        print $O2 "$chr_seg\t$random_start\t$random_end\n";
                        # print "$.\t$chr_seg\t$start\t$end\n";
                        # print "$chr_seg\t$random_start\t$random_end\n";
                    }
                    print "finish\t$i\t$Chr\t$tissue\t$Cutoff\n";
                }
            }
            close($O2);  
        # }
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