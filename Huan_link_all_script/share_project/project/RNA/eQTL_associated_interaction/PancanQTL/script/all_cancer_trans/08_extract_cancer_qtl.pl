#对"../../data/trans_eQTLs_all_re.gz"进行整理，得cancer 对应的文件""../../output/${cancer}/Trans_eQTL/hotspot_trans_eQTL/interval_${j}/qtl_gene.bed.gz";,得hotspot_gene文件${out_dir}/hotspot_gene.bed.gz
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use Parallel::ForkManager;

my $f1 = "../../data/trans_eQTLs_all_re.gz";
# open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件 
my %hash1;
while(<$I1>)
{
    chomp;
    unless(/^cancer_type/){
       my @f= split/\t/;
       my $cancer=$f[0];
       my $chr = $f[2];
       my $pos = $f[3];
       my $start = $pos;
       my $end = $start +1;
       my $gene = $f[5];
       my $p_value = $f[-1];
       if($p_value <5e-8){
            my $k = $cancer;
            my $v= "$chr\t$start\t$end\t$gene\t$p_value\t$cancer";
            push @{$hash1{$k}},$v;
       }

    }
}    

my $j = 18;
my $cutoff = 0.176;
# my @cancers =  ("ACC","BRCA","COAD","ESCA","KICH","KIRC","KIRP","LAML","LIHC","LUAD","LUSC","OV","PAAD","PRAD","SKCM","STAD","TGCT","THCA","UCEC","UCS"); 
my @cancers =  ("BRCA","KIRP"); 
foreach my $k(@cancers){
    my $out_dir = "../../output/${k}/Trans_eQTL/hotspot_trans_eQTL/interval_${j}";
    if(exists $hash1{$k}){
        my @vs = @{$hash1{$k}};
        my $out = join("\n",@vs);
        my $fo1 = "$out_dir/qtl_gene.bed.gz";
        open my $O1, "| gzip >$fo1" or die $!;
        print $O1 "$out\n";
        close($O1);
        my $hotspot = "${out_dir}/${k}_segment_hotspot_cutoff_${cutoff}.bed.gz";
        system "bedtools intersect -a $hotspot -b $fo1  -wa -wb |cut -f1-3,7,9| sort -u |gzip >${out_dir}/hotspot_gene.bed.gz";
    }
}
