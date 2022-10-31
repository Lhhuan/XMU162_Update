#用01_hotspot_egene_idx.txt.gz 为"/home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge_egene_5e_8.bed.gz"添加idx,得"../output/02_hotspot_egene_idx.txt.gz"
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::MoreUtils ':all';

for my $i(1..22){
    my $f1 = "../output/01_chr${i}_hotspot_egene_idx.txt.gz";
    open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
    my $f3 = "/home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge_egene_5e_8.bed.gz";
    open( my $I3 ,"gzip -dc $f3|") or die ("can not open input file '$f3' \n"); #读压缩文件

    my $fo1 = "../output/02_chr${i}_hotspot_egene_idx.txt.gz"; #
    open my $O1, "| gzip >$fo1" or die $!;
    print $O1 "Hotspot_idx\tEgene_idx\n";

    my (%hash1,%hash2,%hash3);

    while(<$I1>)
    {
        chomp;
        unless(/^name/){
            my @f =split/\t/;
            my $name  =$f[0];
            my $idx=$f[1];
            $hash1{$name}=$idx;
        }
    }
    my $chr="chr${i}";
    while(<$I3>)
    {
        chomp;
        unless(/^hotspot_chr/){
            my @f =split/\t/;
            if($f[0]=~/^$chr$/){
                my $hotspot=join("_",@f[0..2]);
                my $egene =$f[3];
                if(exists $hash1{$hotspot}){
                    my $hotspot_idx = $hash1{$hotspot};
                    if(exists $hash1{$egene}){
                        my $egene_idx =$hash1{$egene};
                        print $O1 "$hotspot_idx\t$egene_idx\n";
                    }
                }
            }
        }
    }

    print "$i\n";
}