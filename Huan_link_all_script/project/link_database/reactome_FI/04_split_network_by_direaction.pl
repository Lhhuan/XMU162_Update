#将 ./output/03_FIsInGene_122921_with_annotations_ENSG_filter_predicted.txt 按照方向分为三部分，start-end: "./output/04_ReactomeFI_part1_start_end.txt.gz" end-start: "./output/04_ReactomeFI_part2_end_start.txt.gz" 双向: "./output/04_ReactomeFI_part3_both.txt.gz", 并得无向图"./output/04_ReactomeFI_start_end_all.txt.gz"
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;


my $f1 = "./output/03_FIsInGene_122921_with_annotations_ENSG_filter_predicted.txt";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n";
my $fo1 = "./output/04_ReactomeFI_part1_start_end.txt.gz"; #
# open my $O1, '>', $fo1 or die "$0 : failed to open output file  '$fo1' : $!\n";
open my $O1, "| gzip >$fo1" or die $!;
my $fo2 = "./output/04_ReactomeFI_part2_end_start.txt.gz"; #
open my $O2, "| gzip >$fo2" or die $!;
my $fo3 = "./output/04_ReactomeFI_part3_both.txt.gz"; #
open my $O3, "| gzip >$fo3" or die $!;
my $fo4 = "./output/04_ReactomeFI_start_end_all.txt.gz"; #
open my $O4, "| gzip >$fo4" or die $!;

print $O1 "end\tstart\n";
print $O2 "start\tend\n";
print $O3 "start-end\tstart-end\n";
print $O4 "start\tend\n";



my (%hash1,%hash2);
while(<$I1>)
{
    chomp;
    my @f= split/\t/;
    unless(/^Gene1/){
        my $Direction =$f[3];
        my $Gene1_ensembl =$f[-2];
        my $Gene2_ensembl =$f[-1];
        if ($Direction =~/^<-$/ |$Direction =~/^\|-$/ ){
            print $O2 "$Gene1_ensembl\t$Gene2_ensembl\n";
            print $O4 "$Gene2_ensembl\t$Gene1_ensembl\n";
        }elsif($Direction =~/^->$/ |$Direction =~/^-\|$/){
            
            print $O1 "$Gene1_ensembl\t$Gene2_ensembl\n";
            print $O4 "$Gene1_ensembl\t$Gene2_ensembl\n";
        }
        else{
            print "$Direction\n";
            print $O3 "$Gene1_ensembl\t$Gene2_ensembl\n";
            print $O4 "$Gene1_ensembl\t$Gene2_ensembl\n";   
        }
    }
}
