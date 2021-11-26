# 相同组织，cancer 和 normal -f 1 intersect ，找全部share的hotspot,得my $out1 = "../../output/cancer_total/share/total/${cancer}_contain_${tissue}.bed";my $out2 = "../../output/cancer_total/share/total/${tissue}_contain_${cancer}.bed";
#得汇总文件"../../output/cancer_total/share/total/05_cancer_tissue_intersect_total.bed.gz", cancer 被完全包含的文件"../../output/cancer_total/share/total/05_cancer_share_total.bed.gz"
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

my $fo1 = "../../output/cancer_total/share/total/05_cancer_tissue_intersect_total.bed.gz";
open my $O1, "| gzip >$fo1" or die $!;
my $fo2 = "../../output/cancer_total/share/total/05_cancer_share_total.bed.gz";
open my $O2, "| gzip >$fo2" or die $!;

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
            my $out1 = "../../output/cancer_total/share/total/${cancer}_contain_${tissue}.bed";
            my $out2 = "../../output/cancer_total/share/total/${tissue}_contain_${cancer}.bed";
            my $command1 = "bedtools intersect -f 1 -a $file1 -b $file2 -wo >$out1";
            my $command2 = "bedtools intersect -f 1 -a $file2 -b $file1 -wo >$out2";
            system $command1;
            system $command2;
            
            #_--------------------
            my $f2 = $out1;
            open my $I2, '<', $f2 or die "$0 : failed to open input file '$f2' : $!\n";  

            while(<$I2>)
            {
                chomp;
                my @f= split/\t/;
                print $O1 "$_\t$cancer\t$tissue\n";
                print $O2 "$_\t$cancer\t$tissue\n";
            }   

            print "$cancer\t$tissue\n";

            my $f3 = $out2;
            open my $I3, '<', $f3 or die "$0 : failed to open input file '$f3' : $!\n"; 

            while(<$I3>)
            {
                chomp;
                my @f= split/\t/;
                my $tissue_h = join("\t",@f[0..2]);
                my $cancer_h = join("\t",@f[3..5]);
                my $overlap_bp = $f[-1];
                print $O1 "$cancer_h\t$tissue_h\t$overlap_bp\t$cancer\t$tissue\n";#和 $I2输出相同顺序
                # print $O1 "$_\t$tissue\t$cancer\n";
            }             
            print "$tissue\t$cancer\n";

        }
    }
}





