#利用../../data/pancanQTL_gtex_eQTL.txt，../../output/cancer_total/share/all/${cancer}_contain_${tissue}.bed，"../../output/${cancer}/Cis_eQTL/hotspot_cis_eQTL/interval_${j}/${cancer}_segment_hotspot_cutoff_${cutoff}.bed.gz";"../../../GTEx/output/${tissue}/Cis_eQTL/hotspot_cis_eQTL/interval_${j}/${tissue11}_segment_hotspot_cutoff_${cutoff}.bed.gz" 鉴别组织特异的hotspot"$output_dir/tissue_${cancer}_${tissue}_specific.bed.gz"和cancer特异的hotspot  $output_dir/cancer_${cancer}_${tissue}_specific.bed.gz,cancer 特异汇总文件为../../output/cancer_total/specific/pure/cancer_specific.bed.gz，tissue特异汇总文件../../output/cancer_total/specific/pure/tissue_specific.bed.gz
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
my $f1 = "../../data/pancanQTL_gtex_eQTL.txt";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 

my $fo3 = "../../output/cancer_total/specific/pure/cancer_specific.bed.gz";
open my $O3, "| gzip >$fo3" or die $!;

print $O3 "h_chr\th_start\th_end\tcancer\ttissue\n";

#--------------------
my $fo4 = "../../output/cancer_total/specific/pure/tissue_specific.bed.gz";
open my $O4, "| gzip >$fo4" or die $!;

print $O4 "h_chr\th_start\th_end\tcancer\ttissue\n";

#------------------------

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
            my $out1 = "../../output/cancer_total/share/all/${cancer}_contain_${tissue}.bed";
            
            my $output_dir = "../../output/cancer_total/specific/pure/${cancer}";
            unless(-e $output_dir){
                system "mkdir -p $output_dir";
            }
            #_--------------------

            my $fo1 = "$output_dir/cancer_${cancer}_${tissue}_specific.bed.gz";
            open my $O1, "| gzip >$fo1" or die $!;
            my $fo2 = "$output_dir/tissue_${cancer}_${tissue}_specific.bed.gz";
            open my $O2, "| gzip >$fo2" or die $!;


            my $f2 = $out1;
            open my $I2, '<', $f2 or die "$0 : failed to open input file '$f2' : $!\n";  
            my (%hash2,%hash3);
            while(<$I2>)
            {
                chomp;
                my @f= split/\t/;
                my $C_chr =$f[0];
                my $C_start =$f[1];
                my $C_end = $f[2];
                my $N_chr = $f[3];
                my $N_start =$f[4];
                my $N_end = $f[5];
                my $k2 = "$C_chr\t$C_start\t$C_end";
                my $k3 = "$N_chr\t$N_start\t$N_end";
                $hash2{$k2}=1;
                $hash3{$k3}=1;
            }   
            #------------------------------------
            my $f3 = $file1; #cancer
            open( my $I3 ,"gzip -dc $f3|") or die ("can not open input file '$f3' \n"); #读压缩文件  
            while(<$I3>)
            {
                chomp;
                my @f= split/\t/;
                my $chr =$f[0];
                my $start =$f[1];
                my $end = $f[2];
                my $k2 = "$chr\t$start\t$end";
                unless(exists $hash2{$k2}){
                    print $O1 "$_\n";
                    print $O3 "$_\t$cancer\t$tissue\n";
                }
            } 

            my $f4 = $file2; #normal
            open( my $I4 ,"gzip -dc $f4|") or die ("can not open input file '$f4' \n"); #读压缩文件  
            while(<$I4>)
            {
                chomp;
                my @f= split/\t/;
                my $chr =$f[0];
                my $start =$f[1];
                my $end = $f[2];
                my $k3 = "$chr\t$start\t$end"; 
                unless(exists $hash3{$k3}){
                    print $O2 "$_\t$cancer\t$tissue\n";
                    print $O4 "$_\t$cancer\t$tissue\n";
                }
            } 
            #------------------------------------
            print "$cancer\t$tissue\n";
        }
    }
}

