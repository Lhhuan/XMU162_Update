#将../data/* 下面normalizeda好的文件合在一起

#!/usr/bin/perl
use warnings;
use strict; 
use utf8;

my @datasets = ("Arena-Idb","LncACTdb","miRecords","miRTarBase","RISE","RNAInter","starBase");
foreach my $dataset(@datasets){
    my $path = "../data/${dataset}/normalized";
    my $command = "cat ${path}/* >> ../output/01_all_RNA_associated_RNA-RNA_Protein_and_DNA_interaction.txt";
    system "$command";
}

my $f1 = "../output/01_all_RNA_associated_RNA-RNA_Protein_and_DNA_interaction.txt";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n";
# open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
my $fo1 = "../output/01_all_RNA_associated_RNA-RNA_Protein_and_DNA_interaction.txt.gz";
open my $O1, "| gzip >$fo1" or die $!;

my $header = "Chr1\tStart1\tEnd1\tStrand1\tType1\tinteraction_name1\tGene_ID1\tGene_name1\tEntrez_ID1\tChr2\tStart2\tEnd2\tStrand2\tType2\tinteraction_name2\tGene2_ID2\tGene_name2\tEntrez_ID2\tSub_interaction_type\tInteraction_type\tGenome_version\tTissue/Cell_line\tpancancerNum\tSource";
print $O1 "$header\n";

while(<$I1>) #去掉重复的header
{
    chomp;
    unless(/^Chr1/){
        print $O1 "$_\n";

    }
}

close $O1;
system "rm ../output/01_all_RNA_associated_RNA-RNA_Protein_and_DNA_interaction.txt"; #删掉中间文件，节省空间
    