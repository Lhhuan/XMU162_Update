#cancer两两之间进行 intersect ，找全部share的hotspot,得汇总文件"../../output/cancer_total/share/total/31_absolute_cancer_intersect.bed.gz"
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use Parallel::ForkManager;
my @cutoffs;
my $cutoff =0.176;
# my $tissue = "Lung";
my $j = 18;
my @tissues =  ("ACC","BRCA","COAD","ESCA","KICH","KIRC","KIRP","LAML","LIHC","LUAD","LUSC","OV","PAAD","PRAD","SKCM","STAD","TGCT","THCA","UCEC","UCS");

my $fo1 = "../../output/cancer_total/share/total/31_absolute_cancer_intersect.bed.gz";
open my $O1, "| gzip >$fo1" or die $!;

print $O1 "cancer1_chr\tcancer1_start\tcancer1_end\tcancer2_chr\tcancer2_start\tcancer2_end\toverlap_bp\tcancer1\tcancer2\n";

for my $tissue1(@tissues){
    my $tissue11 = $tissue1; 
    $tissue11 =~ s/Whole_Blood/whole_blood/g;  #Whole_Blood dir name 和 file name不一样，所以引入两个变量
    my $file1 = "../../output/${tissue1}/Cis_eQTL/hotspot_cis_eQTL/interval_${j}/${tissue11}_segment_hotspot_cutoff_${cutoff}.bed.gz";
    for my $tissue2(@tissues){
        my $tissue22 = $tissue2;
        $tissue22 =~ s/Whole_Blood/whole_blood/g;
        unless ($tissue1 eq $tissue2){
            my $file2 = "../../output/${tissue2}/Cis_eQTL/hotspot_cis_eQTL/interval_${j}/${tissue22}_segment_hotspot_cutoff_${cutoff}.bed.gz";
            my $command = "bedtools intersect -f 1 -a $file1 -b $file2 -wo >tmp.bed";
            # print "$command\n";
            system $command;
            my $f1 = "tmp.bed";
            open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 

            while(<$I1>)
            {
                chomp;
                my @f= split/\t/;
                print $O1 "$_\t$tissue1\t$tissue2\n";
            }           
            print "$tissue1\t$tissue2\n";

        }
    }
}
