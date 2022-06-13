#对../output/nodes_annotation/${markers}_Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz 格式进行调整，得../output/nodes_annotation/Adjust_${markers}_Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
use Parallel::ForkManager;
use List::MoreUtils ':all';

my @markers = ("H3K27ac","H3K4me1","H3K4me3","H3K9ac","H3K36me3","H3K27me3","H3K9me3","CTCF","CHROMATIN_Accessibility","TFBS");
foreach my $marker(@markers){
    my $f1 = "../output/nodes_annotation/${marker}_Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz";
    open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
    my $fo1 = "../output/nodes_annotation/Adjust_${marker}_Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz"; #
    open my $O1, "| gzip >$fo1" or die $!;
     
    my %hash1;
    print "$marker\n";
    if($marker =~/CHROMATIN_Accessibility/){
        print $O1 "Hotspot\t${marker}\n"; #postion:cell_line
        while(<$I1>)
        {
            chomp;
            my @f =split/\t/;
            my $h_chr=$f[0];
            my $h_start =$f[1];
            my $h_end =$f[2];
            my $f_pos =join("_",@f[3..5]);
            my $cell_line =$f[7];
            my $k= join("_",@f[0..2]);
            my $v = "$f_pos:$cell_line";
            push @{$hash1{$k}},$v;
        }
        foreach my $k(sort keys %hash1){
            my @vs= @{$hash1{$k}};
            @vs = uniq(@vs);
            my $v=join(";",@vs);
            print $O1 "$k\t$v\n";
        }
        close($O1);
    }elsif($marker =~/TFBS/){
        print $O1 "Hotspot\t${marker}\n";#Factor:pos:cell_line
        while(<$I1>)
        {
            chomp;
            my @f =split/\t/;
            my $h_chr=$f[0];
            my $h_start =$f[1];
            my $h_end =$f[2];
            my $f_pos =join("_",@f[3..5]);
            my $factor =$f[6];
            my $cell_line =$f[7];
            my $k= join("_",@f[0..2]);
            my $v="$factor:$f_pos:$cell_line";
            push @{$hash1{$k}},$v;
        }
        foreach my $k(sort keys %hash1){
            my @vs= @{$hash1{$k}};
            @vs = uniq(@vs);
            my $v=join(";",@vs);
            print $O1 "$k\t$v\n";
        }
        close($O1);
    }
    else{
        print $O1 "Hotspot\t${marker}\n";
        while(<$I1>)
        {
            chomp;
            my @f =split/\t/;
            my $h_chr=$f[0];
            my $h_start =$f[1];
            my $h_end =$f[2];
            my $f_pos =join("_",@f[3..5]);
            my $cell_line =$f[6];
            my $k= join("_",@f[0..2]);
            my $v="$f_pos:$cell_line";
            push @{$hash1{$k}},$v;
        }
        foreach my $k(sort keys %hash1){
            my @vs= @{$hash1{$k}};
            @vs = uniq(@vs);
            my $v=join(";",@vs);
            print $O1 "$k\t$v\n";
        }
        close($O1);
    }
}