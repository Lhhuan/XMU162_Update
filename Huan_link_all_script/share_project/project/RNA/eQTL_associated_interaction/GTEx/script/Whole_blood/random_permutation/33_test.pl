# 对 resemble hotspot,"$output_dir/${i}_resemble_${input_file_base_name} 和../../../output/${tissue}/Cis_eQTL/${group}_cis_eQTL/interval_18/whole_blood_segment_${group}_cutoff_${cutoff}.bed.gz 进行 random_permutation，得../../../output/${tissue}/Cis_eQTL/${group}_cis_eQTL/interval_18_random/random_permutation/${cutoff}/${i}_resemble_${tissue}_segment_${group}_cutoff_${cutoff}.bed.gz
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
use Parallel::ForkManager;

my $cutoff = 0.176;
my $group = "hotspot";
my $tissue = "Whole_Blood";


my $f1 = "../../../output/${tissue}/Cis_eQTL/${group}_cis_eQTL/interval_18/whole_blood_segment_${group}_cutoff_${cutoff}.bed.gz"; #hotspot

open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件

my %hash1;
my @hotspot=();
while(<$I1>)
{
    chomp;
    push @hotspot,$_;
}
my $length_hotspot = @hotspot;

for (my $i=1;$i<10;$i++){
    my @hotspot_and_resemble=@hotspot;
#    push @hotspot_and_resemble,@hotspot;
    my $length = @hotspot_and_resemble;
    # print "\n";
    my $f2 = "../../../output/${tissue}/Cis_eQTL/${group}_cis_eQTL/interval_18_random/original_random/${cutoff}/${i}_resemble_${tissue}_segment_${group}_cutoff_${cutoff}.bed.gz"; #resemble
    open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件
    my $fo1 = "../../../output/${tissue}/Cis_eQTL/${group}_cis_eQTL/interval_18_random/random_permutation/${cutoff}/${i}_resemble_${tissue}_segment_${group}_cutoff_${cutoff}.bed.gz";
    open my $O1, "| gzip >$fo1" or die $!;
    while(<$I2>)
    {
        chomp;
        push @hotspot_and_resemble,$_;
    }
    my $num = @hotspot_and_resemble;
    print "$length\t$num\n";
    #-----
    # my @random_permutations = randomElem ( $length_hotspot, @hotspot_and_resemble ) ; # pick any $length_hotspot from @hotspot_and_resemble ，把$num和@array传递给子程序。这里是用的值传递。还有一种方式是引用传递，相当于硬链接
    # foreach my $random_permutation(@random_permutations){
    #     print $O1 "$random_permutation\n";
    # }
    print "$i\n";
}




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