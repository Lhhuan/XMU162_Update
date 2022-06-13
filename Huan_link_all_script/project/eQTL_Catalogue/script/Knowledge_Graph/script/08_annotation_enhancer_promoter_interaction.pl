# 用/home/huanhuan/project/link_database/OncoBase/output/hg38/01_${type}_sorted.bed.gz 对"/home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge_egene_0.05.bed.gz"进行annotation,得"../output/edges_annotation/success_${type}_hotspot_egene.bed.gz"，adjust format得"../output/edges_annotation/${type}_hotspot_egene.bed.gz"
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
use Parallel::ForkManager;
use List::MoreUtils ':all';

my %hash1;
my $sorted_input_file = "/home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge_egene_0.05.bed.gz";


my @types=("TSS_TSS","ENH_ENH","TSS_ENH");
foreach my $type(@types){
    my $fo1 = "../output/edges_annotation/success_${type}_hotspot_egene.bed.gz";
    open my $O1, "| gzip >$fo1" or die $!;
    print $O1 "Hotspot_gene\t${type}\n";
    #-------------
    my $fo2 = "../output/edges_annotation/success_${type}_ori.txt.gz";
    open my $O2, "| gzip >$fo2" or die $!;
    #-------------
    my $f1= "../output/edges_annotation/${type}_hotspot_egene.bed.gz";
    system "bedtools intersect -a $sorted_input_file -b /home/huanhuan/project/link_database/OncoBase/output/hg38/01_${type}_sorted.bed.gz -wa -wb |gzip > ${f1}" ;
    # print $O1 "hchr_hstart_hend\tegene\ttissue:p_value";
    open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
    my %hash2;
    while(<$I1>)
    {
        chomp;
        my @f=split/\t/;
        my $hotspot=join("_",@f[0..2]);
        my $ensg= $f[3];
        my $f_pos =join("_",@f[4..6]);
        my $enhancer_target_gene= $f[-2];
        my $ori_bed =$f[-1];
        my @tt =split/\;/,$ori_bed;
        my $tissue = $tt[6];
        my $tissue_id =$tissue;
        $tissue_id =~s/_.*//g;
        unless($enhancer_target_gene =~/-/){
            my @t =split/:/,$enhancer_target_gene;
            my $en_ensg =$t[0];
            $en_ensg =~ s/\..*//g;
            # print "$t[0]\t$en_ensg\n";
            my $symbol=$t[1];
            if($ensg eq $en_ensg){
                my $k = "$hotspot:$ensg";
                push @{$hash1{$k}},$tissue_id;
                unless(exists $hash2{$ori_bed}){
                    $hash2{$ori_bed}=1;
                    print $O2 "$ori_bed\n";
                }
            }
        }
    }
    foreach my $k(sort keys %hash1){
        my @vs=@{$hash1{$k}};
        @vs=uniq(@vs);
        my $v= join(";",@vs);
        print $O1 "$k\t$v\n";
    }
}

