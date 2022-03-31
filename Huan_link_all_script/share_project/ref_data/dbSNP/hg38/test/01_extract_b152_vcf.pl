#将GCF_000001405.25.bgz 中所需要的 id,chr,pos,ref,alt 提取出来，得01_extract_b152_vcf.txt.gz
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;

my $f1 = "GCF_000001405.25.bgz";
# open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n";
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
my $fo1 = "01_extract_b152_vcf.txt.gz";
open my $O1, "| gzip >$fo1" or die $!;

print $O1 "Chr\tPos\tRS_id\tRef\tAlt\n";
while(<$I1>) #去掉重复的header
{
    chomp;
    my @f = split/\s+/;
    unless(/^#/){
        my $chr = $f[0];
        $chr =~ s/NC_000001\.//g;
        # print "$chr\n";
        my $pos =$f[1];
        my $Rs_id = $f[2];
        my $REF =$f[3];
        my $ALT =$f[4];
        my $output = "$chr\t$pos\t$Rs_id\t$REF\t$ALT";
        print $O1 "$output\n";
    }
}

close $O1;