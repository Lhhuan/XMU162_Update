# 相同组织，cancer 和 normal  intersect ，找全部overlap(为找cancer specifc做铺垫)的hotspot,得my $out1 = "../../output/cancer_total/share/all/${cancer}_contain_${tissue}.bed";
#得汇总文件"../../output/cancer_total/share/all/06_cancer_tissue_intersect_all.bed.gz",
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

my $out_dir = "../../output/cancer_total/share/all";

my $fo1 = "$out_dir/06_cancer_tissue_intersect_all.bed.gz";
open my $O1, "| gzip >$fo1" or die $!;


my $f1 = "../../data/pancanQTL_gtex_eQTL.txt";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
my %hash1;
while(<$I1>)
{
    chomp;
    unless(/^Study/){
       my @f= split/\t/;
       my $cancer=$f[0];
       $cancer =~ s/\s+//g;
       my $GTEx_tissue=$f[-1];
       my @ts=split/;/,$GTEx_tissue;
    #    print "$cancer\n";
       foreach my $t(@ts){
           push @{$hash1{$cancer}},$t;
       }
    }
}    

my @cancers =  ("ACC","BRCA","COAD","ESCA","KICH","KIRC","KIRP","LAML","LIHC","LUAD","LUSC","OV","PAAD","PRAD","SKCM","STAD","TGCT","THCA","UCEC","UCS"); 

foreach my $cancer(@cancers){
    # print "$cancer\n";
    if(exists $hash1{$cancer}){
        # print "$cancer\n";
        my @tissues =@{$hash1{$cancer}};
        foreach my $tissue(@tissues){
            my $file1 = "../../output/${cancer}/Cis_eQTL/hotspot_cis_eQTL/interval_${j}/${cancer}_segment_hotspot_cutoff_${cutoff}.bed.gz";
            my $tissue11 = $tissue;
            $tissue11 =~ s/Whole_Blood/whole_blood/g;  #Whole_Blood dir name 和 file name不一样，所以引入两个变量
            my $file2 = "../../../GTEx/output/${tissue}/Cis_eQTL/hotspot_cis_eQTL/interval_${j}/${tissue11}_segment_hotspot_cutoff_${cutoff}.bed.gz";
            my $out1 = "${out_dir}/${cancer}_contain_${tissue}.bed";
            # my $out2 = "${out_dir}/${tissue}_contain_${cancer}.bed";
            my $command1 = "bedtools intersect  -a $file1 -b $file2 -wo > $out1";
            
            system $command1;
            #_--------------------
            my $f2 = $out1;
            open my $I2, '<', $f2 or die "$0 : failed to open input file '$f2' : $!\n";  

            while(<$I2>)
            {
                chomp;
                my @f= split/\t/;
                print $O1 "$_\t$cancer\t$tissue\n";
            }   

            print "$cancer\t$tissue\n";
        }
    }
}

