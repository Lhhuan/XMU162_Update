#对../../data/cis_eQTLs_all_re.gz进行整理，得cancer 对应的文件"../../output/cancer_total/specific/pure/${cancer}/qtl_gene.bed.gz";
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use Parallel::ForkManager;

my $f1 = "../../data/cis_eQTLs_all_re.gz";
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
    #    if($p_value <5e-8){
            my $k = $cancer;
            my $v= "$chr\t$start\t$end\t$gene\t$p_value\t$cancer";
            push @{$hash1{$k}},$v;
    #    }

    }
}    

my @cancers =  ("ACC","BRCA","COAD","ESCA","KICH","KIRC","KIRP","LAML","LIHC","LUAD","LUSC","OV","PAAD","PRAD","SKCM","STAD","TGCT","THCA","UCEC","UCS"); 

foreach my $k(@cancers){
    if(exists $hash1{$k}){
        my @vs = @{$hash1{$k}};
        my $out = join("\n",@vs);
        my $fo1 = "../../output/cancer_total/specific/pure/${k}/qtl_gene.bed.gz";
        open my $O1, "| gzip >$fo1" or die $!;
        print $O1 "$out\n";
    }
}

