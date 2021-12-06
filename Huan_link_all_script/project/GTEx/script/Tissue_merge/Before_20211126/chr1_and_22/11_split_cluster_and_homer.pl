
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
use Parallel::ForkManager;

# my %hash1;

for(my $i=0;$i<6;$i++){
    # print "$i\n";
    my $f1 = "/share/data0/QTLbase/huan/GTEx/Tissue_merge/Cis_eQTL/hotspot/interval_18/kmer/6/chr1_and_22/communities.bed.gz";
    # open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
    open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
    # my $fo1 = "/share/data0/QTLbase/huan/GTEx/Tissue_merge/Cis_eQTL/hotspot/interval_18/kmer/6/chr1_and_22/homer/communities_${i}.bed.gz";
    # open my $O1, "| gzip >$fo1" or die $!;  
    my $fo1 = "/share/data0/QTLbase/huan/GTEx/Tissue_merge/Cis_eQTL/hotspot/interval_18/kmer/6/chr1_and_22/homer/communities_${i}.bed";
    open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
    while(<$I1>)
    {
        chomp;
        my @f= split/\t/;
        my $chr = $f[0];
        my $start = $f[1];
        my $end =$f[2];
        my $cluster = $f[3];
        # print "$cluster\n";
        if($cluster =~/$i/){
            my $output = join("\t",@f[0..2]);
            print $O1 "$output\n";
        }
    } 
    close($O1);
    my $output_dir = "/share/data0/QTLbase/huan/GTEx/Tissue_merge/Cis_eQTL/hotspot/interval_18/kmer/6/chr1_and_22/homer/${i}/";
    if(-e $output_dir){
        print "${output_dir}\texist\n";
    }
    else{
        system "mkdir -p $output_dir";
    }

    system "findMotifsGenome.pl $fo1 /home/huanhuan/ref_data/gencode/GRCh37.primary_assembly.genome.fa $output_dir -size 200";
    print "$i finish\n";
}

# findMotifsGenome.pl /share/data0/QTLbase/huan/GTEx/Tissue_merge/Cis_eQTL/hotspot/interval_18/kmer/6/chr1_and_22/homer/communities_0.bed ~/ref_data/gencode/GRCh37.primary_assembly.genome.fa /share/data0/QTLbase/huan/GTEx/Tissue_merge/Cis_eQTL/hotspot/interval_18/kmer/6/chr1_and_22/homer/0/ -size 200