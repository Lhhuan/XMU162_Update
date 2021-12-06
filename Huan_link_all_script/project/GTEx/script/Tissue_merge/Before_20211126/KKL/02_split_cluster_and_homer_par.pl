
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
use Parallel::ForkManager;

# my %hash1;
my $pm = Parallel::ForkManager->new(10); 
for(my $i=1;$i<7;$i++){
# for(my $i=1;$i<2;$i++){
    my $pid = $pm->start and next; #开始多线程
    # print "$i\n";
    my $f1 = "/home/huanhuan/project/GTEx/script/Tissue_merge/KKL/01_kkl_result.txt";
    open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
    # open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
    # my $fo1 = "/share/data0/QTLbase/huan/GTEx/Tissue_merge/Cis_eQTL/hotspot/interval_18/kmer/6/chr1_and_22/homer/communities_${i}.bed.gz";
    # open my $O1, "| gzip >$fo1" or die $!;  
    my $fo1 = "/home/huanhuan/project/GTEx/script/Tissue_merge/KKL/homer/communities_${i}.bed";
    open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
    while(<$I1>)
    {
        chomp;
        unless(/^name/){
            my @f= split/\t/;
            my $hotspot=$f[0];
            my $cluster = $f[1];
            my @t=split/\:/,$hotspot;
            my $chr= $t[0];
            my @s=split/\-/,$t[1];
            my $start =$s[0];
            my $end=$s[1];
            # print "$cluster\n";
            if($cluster =~/$i/){
                print $O1 "$chr\t$start\t$end\n";
            }
        }
    } 
    close($O1);
    my $output_dir = "/home/huanhuan/project/GTEx/script/Tissue_merge/KKL/homer/${i}/";
    if(-e $output_dir){
        print "${output_dir}\texist\n";
    }
    else{
        system "mkdir -p $output_dir";
    }

    system "findMotifsGenome.pl $fo1 hg19 $output_dir -size given";
    print "$i finish\n";
    $pm->finish;  #多线程结束
}

# findMotifsGenome.pl /share/data0/QTLbase/huan/GTEx/Tissue_merge/Cis_eQTL/hotspot/interval_18/kmer/6/chr1_and_22/homer/communities_0.bed ~/ref_data/gencode/GRCh37.primary_assembly.genome.fa /share/data0/QTLbase/huan/GTEx/Tissue_merge/Cis_eQTL/hotspot/interval_18/kmer/6/chr1_and_22/homer/0/ -size 200