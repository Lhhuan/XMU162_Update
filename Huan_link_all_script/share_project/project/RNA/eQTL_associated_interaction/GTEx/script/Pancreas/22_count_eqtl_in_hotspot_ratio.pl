#统计/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/${tissue}/Cis_eQTL/hotspot_cis_eQTL_overlap_SNP/${tissue}_segment_hotspot_cutoff_${cutoff}_SNP.bed.gz中eQTL在SNP中的比例得/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/${tissue}/Cis_eQTL/22_${tissue}_cis_eQTL_in_hotspot_ratio.gz
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use Parallel::ForkManager;


my $tissue= "Liver";
my @cutoffs = ();
for (my $i=0.05;$i<0.31;$i=$i+0.01){ #对文件进行处理，把所有未定义的空格等都替换成NONE
    push @cutoffs,$i;
}


my $fo1 = "/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/${tissue}/Cis_eQTL/22_${tissue}_cis_eQTL_in_hotspot_ratio.gz";
open my $O1, "| gzip >$fo1" or die $!;

print $O1 "Cutoff\tHotspot_chr\tHotspot_start\tHotspot_end\tAll_eQTL_in_hotspot_count\tAll_SNP_in_hotspot_count\teQTL_ratio_in_hotspot\n";

foreach my $cutoff(@cutoffs){
    my (%hash1,%hash2);
    my $input_dir = "/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/${tissue}/Cis_eQTL/hotspot_cis_eQTL_overlap_SNP";
    my $f1 = "${input_dir}/${tissue}_segment_hotspot_cutoff_${cutoff}_SNP.bed.gz";
    # open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
    open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件

    while(<$I1>)
    {
        chomp;
        my @f = split/\t/;
        my $chr = $f[0];
        my $SNP_pos = $f[1];
        my $Pvalue =$f[3];
        my $hotspot_start = $f[-2];
        my $hotspot_end = $f[-1];
        my $k_hotspot = "$chr\t$hotspot_start\t$hotspot_end";
        my $snp_v =  "$chr\t$SNP_pos";
        push @{$hash1{$k_hotspot}},$snp_v;
        if ($Pvalue <5e-8){ #----significant eQTL
            push @{$hash2{$k_hotspot}},$snp_v;
        }
    }
    
    foreach my $k (sort keys %hash1){
        my @all_snp = @{$hash1{$k}};
        my %hash3;
        @all_snp = grep { ++$hash3{$_} < 2 } @all_snp;
        my $all_snp_in_hotspot_num = @all_snp;
        my @eQTL= ();
        if (exists $hash2{$k}){
            my @eqtl = @{$hash2{$k}};
            my %hash4;
            @eqtl = grep { ++$hash4{$_} < 2 } @eqtl;
            my $all_eqtl_in_hotspot_num = @eqtl;
            push @eQTL,$all_eqtl_in_hotspot_num;
        }
        else{
            my $all_eqtl_in_hotspot_num=0;
            push @eQTL,$all_eqtl_in_hotspot_num;
        }
        #--------------------
        my $all_eqtl_in_hotspot_num = $eQTL[0];
        my $ratio = $all_eqtl_in_hotspot_num/$all_snp_in_hotspot_num;
        print $O1 "$cutoff\t$k\t$all_eqtl_in_hotspot_num\t$all_snp_in_hotspot_num\t$ratio\n";
    }
    print "$cutoff\n";
    close($I1);
}

close($O1);
