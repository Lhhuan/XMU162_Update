#adjust ../output/test/021_hotsopt_valid.csv.gz to json得"../output/test/02_valid_hotspot.json.gz"，"../output/02_all_hotspot.csv.gz"减去"../output/test/021_hotsopt_valid.csv.gz"得 "../output/train/02_train_hotspot.csv.gz"
#!/usr/bin/perl
use warnings;
use strict; 
# use utf8;
# binmode STDOUT, ":utf8";
use JSON;

my $f1 = "../output/test/021_chr1_and_22_hotsopt_valid.csv.gz";
# open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n";
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
my $f2 = "../output/02_chr1_22_hotspot.csv.gz";
open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件
my $fo1 = "../output/test/02_valid_hotspot.json"; #
open my $O1, '>', $fo1 or die "$0 : failed to open output file  '$fo1' : $!\n";
# open my $O1, "| gzip >$fo1" or die $!;
# my $fo2 = "../output/train/02_train_hotspot.csv.gz"; #
# open my $O2, "| gzip >$fo2" or die $!;
my $fo2 = "../output/train/02_train_hotspot.csv"; #
open my $O2, '>', $fo2 or die "$0 : failed to open output file  '$fo2' : $!\n";
# open my $O1, '>', $fo1 or die "$0 : failed to open output file  '$fo1' : $!\n";
print $O2 ",Source node name,Target node name,Source node type,Target node type,Relationship type\n";

my (%hash1,%hash2);
while(<$I1>)
{
    chomp;
    unless(/^Source/){
        my @f= split/\,/;
        $hash1{$_}=1;
        my $gene =$f[0];
        my $hotspot= $f[1];
        # print "$gene\t$hotspot\n";
        push @{$hash2{$hotspot}},$gene;
    }
}


my $json = encode_json \%hash2;
print $O1 $json;

while(<$I2>)
{
    chomp;
    unless(/^Source/){
        my @f= split/\,/;
        my $gene =$f[0];
        unless(exists $hash1{$_}){
            print $O2 "$.,$_\n";
        }
    }
}




