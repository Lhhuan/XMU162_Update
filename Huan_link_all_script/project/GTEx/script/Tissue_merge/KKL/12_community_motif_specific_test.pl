#统计每个community中各个部分显著motif"${output_dir}/knowResult_merge.txt"， overlap情况得"${output_dir}/knowResult_merge_class_merge.txt",overlap ratio 文件"./6kmer/community_motif_overlap_ratio.txt"
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
use Parallel::ForkManager;

# my %hash1;
# my @aa=(3..6);

# my $fo3 = "./6kmer/community_motif_overlap_ratio.txt";
# open my $O3, '>', $fo3 or die "$0 : failed to open output file '$fo3' : $!\n";
# print $O3 "Number_of_all_motif\tNumber_of_overlap_motif\toverlap_ratio\tcommunity\n";

# foreach my $j(@aa){

    my $output_dir = "/home/huanhuan/project/GTEx/script/Tissue_merge/chr1_6/6kmer/5_community/";
    if(-e $output_dir){
        print "$output_dir\texist\n";
    }
    else{
        system "mkdir -p $output_dir";
    }
    my $fo1 = "${output_dir}/knowResult_merge1.txt";
    open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
    my $fo2 = "${output_dir}/knowResult_specific.txt";
    open my $O2, '>', $fo2 or die "$0 : failed to open output file '$fo2' : $!\n";
    print $O1 "Motif_Name\tConsensus\tq_value\tcommunity\n";
    print $O2 "Motif_Name\tConsensus\tcommunity\n";
    for(my $i=1;$i<3;$i++){
        # print "$i\n";
        my $f1 = "/share/data0/QTLbase/huan/GTEx/Tissue_merge/Cis_eQTL/hotspot/interval_18/chr1_6/kmer/6/3_community/homer/${i}/knownResults.txt";
        open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
        # open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
        # my $fo1 = "/share/data0/QTLbase/huan/GTEx/Tissue_merge/Cis_eQTL/hotspot/interval_18/kmer/6/chr1_and_22/homer/communities_${i}.bed.gz";
        # open my $O1, "| gzip >$fo1" or die $!;  

        while(<$I1>)
        {
            chomp;
            unless(/^Motif/){
                my @f= split/\t/;    
                my $motif =$f[0] ;
                my $Consensus =$f[1];
                my $q_value = $f[4];
                if ($q_value <0.05){
                    print $O1 "$motif\t$Consensus\t$q_value\t$i\n";
                    # my $v=  
                }
            }

        } 
        print "$i\n";
    }
    for(my $i=3;$i<4;$i++){
        # print "$i\n";
        my $f1 = "/share/data0/QTLbase/huan/GTEx/Tissue_merge/Cis_eQTL/hotspot/interval_18/chr1_6/kmer/6/4_community/homer/${i}/knownResults.txt";
        open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
        # open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
        # my $fo1 = "/share/data0/QTLbase/huan/GTEx/Tissue_merge/Cis_eQTL/hotspot/interval_18/kmer/6/chr1_and_22/homer/communities_${i}.bed.gz";
        # open my $O1, "| gzip >$fo1" or die $!;  

        while(<$I1>)
        {
            chomp;
            unless(/^Motif/){
                my @f= split/\t/;    
                my $motif =$f[0] ;
                my $Consensus =$f[1];
                my $q_value = $f[4];
                if ($q_value <0.05){
                    print $O1 "$motif\t$Consensus\t$q_value\t$i\n";
                    # my $v=  
                }
            }

        } 
        print "$i\n";
    }
    for(my $i=4;$i<6;$i++){
        # print "$i\n";
        my $f1 = "/share/data0/QTLbase/huan/GTEx/Tissue_merge/Cis_eQTL/hotspot/interval_18/chr1_6/kmer/6/5_community/homer/${i}/knownResults.txt";
        open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
        # open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
        # my $fo1 = "/share/data0/QTLbase/huan/GTEx/Tissue_merge/Cis_eQTL/hotspot/interval_18/kmer/6/chr1_and_22/homer/communities_${i}.bed.gz";
        # open my $O1, "| gzip >$fo1" or die $!;  

        while(<$I1>)
        {
            chomp;
            unless(/^Motif/){
                my @f= split/\t/;    
                my $motif =$f[0] ;
                my $Consensus =$f[1];
                my $q_value = $f[4];
                if ($q_value <0.05){
                    print $O1 "$motif\t$Consensus\t$q_value\t$i\n";
                    # my $v=  
                }
            }

        } 
        print "$i\n";
    }

    close($O1);


    my $f2 = $fo1;
    open my $I2, '<', $f2 or die "$0 : failed to open input file '$f2' : $!\n"; 

    my %hash1;
    while(<$I2>)
    {
        chomp;
        unless(/^Motif/){
            my @f= split/\t/;    
            my $motif =$f[0] ;
            my $Consensus =$f[1];
            my $q_value = $f[2];
            my $class = $f[-1];
            # my $k="$Consensus";
            # my $k="$motif";
            my $k="$motif\t$Consensus";
            push @{$hash1{$k}},$class;
        }
    } 
    my %hash3;
    foreach my $k(sort keys %hash1){
        my @vs = @{$hash1{$k}};
        my %hash2;
        @vs = grep { ++$hash2{$_} < 2 } @vs;
        my $number = @vs;
        if($number <2){
            print $O2 "$k\t$vs[0]\n";
            print "@vs\n";
        }
    }
    close($O2);
# }









# findMotifsGenome.pl /share/data0/QTLbase/huan/GTEx/Tissue_merge/Cis_eQTL/hotspot/interval_18/kmer/6/chr1_and_22/homer/communities_0.bed ~/ref_data/gencode/GRCh37.primary_assembly.genome.fa /share/data0/QTLbase/huan/GTEx/Tissue_merge/Cis_eQTL/hotspot/interval_18/kmer/6/chr1_and_22/homer/0/ -size 200