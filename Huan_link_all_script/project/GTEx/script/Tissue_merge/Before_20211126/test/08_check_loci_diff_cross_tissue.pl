#将每个组织进行显著eQTL-snp的chr1进行提取，得"./tmp_output/07_chr_muliti_tissue_eqtl.bed.gz， muti-tissue hotspot的chr1进行提取得 ./tmp_output/07_chr_hotspot.bed.gz
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;
use Parallel::ForkManager;


my $f1 = "./tmp_output/07_chr_muliti_tissue_eqtl_hotspot.bed.gz";
# open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
my $fo1 = "./tmp_output/08_chr_hotspot_ref_alt_uniformity.bed.gz";
# open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
open my $O1, "| gzip >$fo1" or die $!;
my $fo2 = "./tmp_output/08_chr_hotspot_ref_alt_diff.bed.gz";
# # open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
open my $O2, "| gzip >$fo2" or die $!;
# print "$tissue\tstart\n";
my (%hash1,%hash2);
while(<$I1>)
{
    chomp;
    my @f =split/\t/;
    my $hotspot = join("\t",@f[0..2]);
    my $qtl_chr =  $f[3];
    my $qtl_start =  $f[4];
    my $qtl_end =  $f[5];  
    my $ref = $f[6];
    my $alt = $f[7];
    my $tissue = $f[8];
    my $k1 = join("\t",@f[0..5]);
    my $v1 = "$ref\t$alt";
    my $k2 = "$k1\t$v1";
    my $v2 = $tissue;
    push @{$hash1{$k1}},$v1;
    push @{$hash2{$k2}},$v2;
}

foreach my $k1(sort keys %hash1){
    my @v1s = @{$hash1{$k1}};
    my (%hash3,%hash4);
    @v1s = grep { ++$hash3{$_} < 2 } @v1s;
    my $length= @v1s;
    if($length <2){
        print $O1 "$k1\t$v1s[0]\n";
    }
    else{
        print "$length\n";
        foreach my $v1(@v1s){
            my $k2 = "$k1\t$v1";
            my @v2s =@{$hash2{$k2}};
            foreach my $v2(@v2s){
                print $O2 "$k2\t$v2\n";
            }
        }
    }
}