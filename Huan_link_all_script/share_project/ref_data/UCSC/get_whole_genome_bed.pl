#过滤hg19.chrom.sizes中chr1_22得hg19.chrom1_22.sizes
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;


my $fo1 = "hg19.chrom1_22_sizes_sorted.bed";
open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
# open my $O1, "| gzip >$fo1" or die $!;

my $f1 = "hg19.chrom1_22_sizes_sorted.txt";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
# open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件

my %hash1;
my %hash2;

while(<$I1>)
{
    chomp;
    my @f = split/\t/;
    my $chr = $f[0];
    my $end = $f[1];
    print $O1 "$chr\t0\t$end\n";
    
}

# foreach my $k(sort keys %hash1){
#     if( exists($hash2{$k})){
#         my $v = $hash2{$k};
#         print $O1 "$v\n";
#     }
# # }



# for (my $i=1;$i<23;$i++){
#     my $k = "chr${i}";
#     $hash1{$k}=1;
# }


# while(<$I1>)
# {
#     chomp;
#     my @f = split/\t/;
#     my $chr = $f[0];
#     $hash2{$chr}=$_;
# }

# foreach my $k(sort keys %hash1){
#     if( exists($hash2{$k})){
#         my $v = $hash2{$k};
#         print $O1 "$v\n";
#     }
# }