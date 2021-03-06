#得100%绝对overlap文件"../../output/cancer_total/21_cancer_pair_ovelap_absolute.txt.gz"，得没有在"../../output/cancer_total/21_cancer_pair_ovelap_absolute.txt.gz"中的片段"../../output/cancer_total/21_cancer_pair_out_ovelap_absolute.txt.gz"
#得用于计算类似jaccard index 数据"../../output/cancer_total/21_cancer_pair_overlap_index_and_related_number.txt.gz"，得类似jaccard index "../../output/cancer_total/21_cancer_pair_overlap_index.txt.gz"
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



my $fo1 = "../../output/cancer_total/21_cancer_pair_ovelap_absolute.txt.gz";
open my $O1, "| gzip >$fo1" or die $!;
print $O1 "Chr1\tstart1\tend1\tChr2\tstart2\tend2\tcancer1\tcancer2\n";

my $fo2 = "../../output/cancer_total/21_cancer_pair_out_ovelap_absolute.txt.gz";
open my $O2, "| gzip >$fo2" or die $!;
print $O2 "cancer1_Chr\tcancer1_start\tcancer1_end\tcancer1\tcancer2\n";

my $fo3 = "../../output/cancer_total/21_cancer_pair_overlap_index_and_related_number.txt.gz";
open my $O3, "| gzip >$fo3" or die $!;
print $O3 "Index\tnumber_of_overlap\tnumber_of_cancer1_out_absolute\tnumber_of_cancer2_out_absolute\tcancer1\tcancer2\n";

my $fo4 = "../../output/cancer_total/21_cancer_pair_overlap_index.txt.gz";
open my $O4, "| gzip >$fo4" or die $!;
print $O4 "Index\tcancer1\tcancer2\n";
my (%hash1,%hash4);
foreach my $tissue1(@tissues){
    # my $file1 = "../../output/${cancer}/Cis_eQTL/hotspot_cis_eQTL/interval_${j}/${cancer}_segment_hotspot_cutoff_${cutoff}.bed.gz";
    my $file1 = "../../output/${tissue1}/Cis_eQTL/hotspot_cis_eQTL/interval_${j}/${tissue1}_segment_hotspot_cutoff_${cutoff}.bed.gz";
    foreach my $tissue2(@tissues){
        print "$tissue1\t$tissue2\n";
        my $final_index = "1\t$tissue1\t$tissue1"; # same tissue
        unless(exists $hash4{$final_index}){
            $hash4{$final_index}=1;
            print $O4 "$final_index\n";
        }
        unless ($tissue1 eq $tissue2){
            my $pair= "$tissue1\t$tissue2";
            my $reverse_pair = "$tissue2\t$tissue1";
            if(exists $hash1{$pair} || $hash1{$reverse_pair}){ #防止重复
                print "$pair\texists\n";
            }
            else{
                $hash1{$pair}=1;
                $hash1{$reverse_pair}=1;
                # my $file1 = "../../output/${cancer}/Cis_eQTL/hotspot_cis_eQTL/interval_${j}/${cancer}_segment_hotspot_cutoff_${cutoff}.bed.gz";
                my $file2 = "../../output/${tissue2}/Cis_eQTL/hotspot_cis_eQTL/interval_${j}/${tissue2}_segment_hotspot_cutoff_${cutoff}.bed.gz";
                my $command1 = "bedtools intersect -f 1 -a $file1 -b $file2 -wo >tmp1.bed";
                my $command2 = "bedtools intersect -f 1 -a $file2 -b $file1 -wo >tmp2.bed";
                system $command1;
                system $command2;
                #-------------------------------------------------- find all tissue1 and tissue2 overlap
                my $f1 = "tmp1.bed"; 
                # open( my $I3 ,"gzip -dc $f3|") or die ("can not open input file '$f3' \n"); #读压缩文件  
                open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n";
                #-------------
                my(%hash2,%hash3,%hash5);
                while(<$I1>)
                {
                    chomp;
                    my @f= split/\t/;
                    my $h1 = join("\t",@f[0..2]);
                    my $h2 = join("\t",@f[3..5]);
                    my $k = "$h1\t$h2";
                    $hash2{$k}=1;
                    $hash3{$h1}=1; #---tissue1
                    $hash5{$h2}=1; #---tissue2
                    print $O1 "$k\t$tissue1\t$tissue2\n";
                }   

                my $f2 = "tmp2.bed"; 
                # open( my $I3 ,"gzip -dc $f3|") or die ("can not open input file '$f3' \n"); #读压缩文件  
                open my $I2, '<', $f1 or die "$0 : failed to open input file '$f2' : $!\n";
                #-------------
                while(<$I2>)
                {
                    chomp;
                    my @f= split/\t/;
                    my $h1 = join("\t",@f[0..2]);
                    my $h2 = join("\t",@f[3..5]);
                    my $k = "$h2\t$h1"; #与tmp1 相反，得到仍然得到tissue1，tissue2;
                    unless(exists $hash2{$k}){
                        $hash2{$k}=1;
                        $hash3{$h2}=1;#---tissue1
                        $hash5{$h1}=1;#---tissue2
                        print $O1 "$k\t$tissue1\t$tissue2\n";
                    }
                } 
                
                #------------------------------------ #-----find tissue1 and tissue2 not in overlap
                my @tissue1_out_absolute=();
                my @tissue2_out_absolute=();
                my $f3 = $file1; #TISSUE1
                open( my $I3 ,"gzip -dc $f3|") or die ("can not open input file '$f3' \n"); #读压缩文件
                while(<$I3>)
                {
                    chomp;
                    my @f= split/\t/;
                    my $h1 = join("\t",@f[0..2]);
                    unless(exists $hash3{$h1}){
                        print $O2 "$h1\t$tissue1\t$tissue2\n";
                        push @tissue1_out_absolute,$h1;
                    }
                } 

                my $f4 = $file2; #TISSUE2
                open( my $I4 ,"gzip -dc $f4|") or die ("can not open input file '$f4' \n"); #读压缩文件
                while(<$I4>)
                {
                    chomp;
                    my @f= split/\t/;
                    my $h1 = join("\t",@f[0..2]);
                    unless(exists $hash5{$h1}){
                        print $O2 "$h1\t$tissue2\t$tissue1\n";
                        push @tissue2_out_absolute,$h1;
                    }
                }

                my $number_of_tissue1_out_absolute = @tissue1_out_absolute;
                my $number_of_tissue2_out_absolute = @tissue2_out_absolute;
                # print $O3 "$number_of_tissue1_out_absolute\t$tissue1\t$number_of_tissue2_out_absolute\t$tissue2\n";
                my $length_of_overlap = keys %hash2;
                my $index = $length_of_overlap/($length_of_overlap+$number_of_tissue1_out_absolute+$number_of_tissue2_out_absolute);
                
                
                print $O3 "$index\t$length_of_overlap\t$number_of_tissue1_out_absolute\t$number_of_tissue2_out_absolute\t$tissue1\t$tissue2\n";
                print $O4 "$index\t$tissue1\t$tissue2\n";


            }
        }
    }
}