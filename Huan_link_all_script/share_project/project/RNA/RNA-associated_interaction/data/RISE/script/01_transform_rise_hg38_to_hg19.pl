#利用../output/RNA1_hg38_to_hg19.bed和../output/RNA2_hg38_to_hg19.bed  将../output/rise_human_all.txt的gene version 转为hg19,得../output/01_rise_human_hg19.txt;
#得因为位置信息不全而不能转出hg19的文件../output/01_rise_human_position_null.txt
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;


my $f1 = "../output/RNA1_hg38_to_hg19.bed";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n";
my $f2 = "../output/RNA2_hg38_to_hg19.bed";
open my $I2, '<', $f2 or die "$0 : failed to open input file '$f2' : $!\n";
my $f3 = "../output/rise_human_all.txt";
open my $I3, '<', $f3 or die "$0 : failed to open input file '$f3' : $!\n";
my $fo1 = "../output/01_rise_human_hg19.txt";
open my $O1, '>', $fo1 or die "$0 : failed to open output file  '$fo1' : $!\n";
my $fo2 = "../output/01_rise_human_position_null.txt"; #得没有匹配的mirna文件
open my $O2, '>', $fo2 or die "$0 : failed to open output file  '$fo2' : $!\n";

# my $header = "Chr1\tStart1\tEnd1\tStrand1\tType1\tProtein1\tGene_ID1\tGene_symbol1\tEntrez_ID1\tChr2\tStart2\tEnd2\tStrand2\tType2\tProtein2\tGene2_ID2\tGene_symbol2\tEntrez_ID2\tSub_interaction_type\tInteraction_type\tGenome_version\tTissue\tpancancerNum\tSource";
# print $O1 "$header\n";
my (%hash1,%hash2,%hash3);

while(<$I1>)
{
    chomp;
    if ($_=~/\->/){
        my @t= split/\->/;
        my $hg38 = $t[0];
        my $hg19 = $t[1];
        $hg19 =~ s/^\s+//g;
        if ($hg19 =~/^chr/){
            $hash1{$hg38}=$hg19;
        }
    }
}
while(<$I2>)
{
    chomp;
    if ($_=~/\->/){
        my @t= split/\->/;
        my $hg38 = $t[0];
        my $hg19 = $t[1];
        $hg19 =~ s/^\s+//g;
        if ($hg19 =~/^chr/){
            $hash2{$hg38}=$hg19;
        }
    }
}

while(<$I3>)
{
    chomp;
    my @f= split/\t/;
    my $number = @f; #22个元素
    if (/^chr1/){
        print $O1 "$_\n";
        print $O2 "$_\n";
    }
    else{
        my $RNA1 = join("\t",@f[1..3]);
        my $RNA2 = join("\t",@f[4..6]);
        my $other_info = join("\t",@f[7..21]);
        my @hg19_RNA1= ();
        my @hg19_RNA2= ();
        if (exists $hash1{$RNA1}){
            my $RNA1_19 = $hash1{$RNA1};
            push @hg19_RNA1, $RNA1_19;
        }
        else{
            print $O2 "$_\n";
        }

        if (exists $hash2{$RNA2}){
            my $RNA2_19 = $hash2{$RNA2};
            push @hg19_RNA2, $RNA2_19;
        }
        else{
            print $O2 "$_\n";
        }
        my $number1= @hg19_RNA1;
        my $number2 = @hg19_RNA2;
        if ($number1>0 && $number2>0){ 
            my $output= "$hg19_RNA1[0]\t$hg19_RNA2[0]\t$other_info";
            print $O1 "$output\n";
        }
    }
}
