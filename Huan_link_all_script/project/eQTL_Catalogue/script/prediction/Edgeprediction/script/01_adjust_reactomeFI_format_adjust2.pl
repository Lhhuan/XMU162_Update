#调整"~/project/link_database/reactome_FI/output/03_FIsInGene_122921_with_annotations_ENSG_filter_predicted.txt"格式得 ../output/train/01_reactome.csv.gz
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;


my $f1 = "/home/huanhuan/project/link_database/reactome_FI/output/03_FIsInGene_122921_with_annotations_ENSG_filter_predicted.txt";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n";
# open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
# my $fo2 = "../output/train/01_reactome.csv.gz"; #
# open my $O2, "| gzip >$fo2" or die $!;

my $fo2 = "../output/train/01_reactome_adjust.csv"; #
open my $O2, '>', $fo2 or die "$0 : failed to open output file  '$fo2' : $!\n";


print $O2 "Source node name,Source node type,Relationship type,Target node type,Target node name\n";

my (%hash1,%hash2);
while(<$I1>)
{
    chomp;
    my @f= split/\t/;
    unless(/^Gene1/){
        my $Gene1 =$f[0];
        my $Gene2 =$f[1];
        my $Gene1_ensembl =$f[-2];
        my $Gene2_ensembl =$f[-1];
        if(exists $hash1{$Gene2} || exists $hash2{$Gene1}){  #判断end是否在start中出现过，或者start 有没有在end出现过，如果出现,交换start和end位置
            print $O2 "$Gene2_ensembl,Gene,ReactomeFI,Gene,$Gene1_ensembl\n";

        }
        else{
            $hash1{$Gene1} =1;
            $hash2{$Gene2}=1;
            print $O2 "$Gene1_ensembl,Gene,ReactomeFI,Gene,$Gene2_ensembl\n";
        }

    }
}

