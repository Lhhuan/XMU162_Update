#分source进行clump,将不同类型QTL按照tissue 和 population 分割成"../output/tissue_pop_used_to_clump/$used_to_clump"，然后用对应的人种文件/share/data0/1kg_phase3_v5_hg19/${Population}/1kg.phase3.v5.shapeit2.${pop}.hg19.all.SNPs.uniq.posID
#进行clump, 得 ../output/tissue_pop_used_to_clump/$used_to_clump.clump

#!/usr/bin/perl
use warnings;
use strict; 
use utf8;


my %hash1;

my $f1 = "/share/data0/QTLbase/data/QTLbase_download_data_sourceid.txt";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n";

while(<$I1>)
{
    chomp;
    unless(/PMID/){
        # print "$file\n";
        my @f =split/\t/;
        my $Sourceid =$f[1];
        my $Tissue =$f[3];
        my $Population =$f[4];
        my $v = "$Tissue\t$Population";
        $hash1{$Sourceid}=$v;
    }
}

#-------------------------------------------------------------#获取组织名称
my $dir = "/share/data0/QTLbase/data/";
opendir (DIR, $dir) or die "can't open the directory!";
my @files = readdir DIR; #获取一个文件夹下的所有文件
my @tissues;
my $suffix = "QTL.txt.gz"; #文件名后缀QTL.txt
foreach my $file(@files){ #分QTL进行
    if ($file =~/$suffix/){ 
        my $f2 = "/share/data0/QTLbase/data/$file";
        open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件

        my %hash2;
        while(<$I2>)
        {
            chomp;
            unless(/SNP_chr/){
                # print "$file\n";
                my @f =split/\t/;
                my $SNP_chr =$f[0];
                my $SNP_pos =$f[1];
                my $Mapped_gene =$f[2];
                my $Trait_chr =$f[3];
                my $Trait_start =$f[4];
                my $Trait_end =$f[5];
                my $Pvalue =$f[-2];
                my $Sourceid =$f[-1];
                if (exists $hash1{$Sourceid}){
                    my $v1 = $hash1{$Sourceid}; # $v1 = "$Tissue\t$Population";
                    my $SNP = "${SNP_chr}_${SNP_pos}";
                    my $P = $Pvalue;
                    my $v2 = "$SNP\t$Mapped_gene\t$Trait_chr\t$Trait_start\t$Trait_end\t$P\t$Sourceid";
                    my $k2 = $v1; # $v1 = "$Tissue\t$Population";
                    push @{$hash2{$k2}},$v2;
                }
            }
        }
        my $type = $file; 
        $type =~ s/\.txt.gz//g;
        foreach my $k2(sort keys %hash2){
            my @vs =  @{$hash2{$k2}};
            my @t = split/\t/,$k2;
            my $Tissue =$t[0];
            my $Population =$t[1];
            my $used_to_clump = "${type}_${Tissue}_${Population}";
            my $fo1 = "../output/tissue_pop_used_to_clump/$used_to_clump";
            open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
            # print $O1 "SNP\tP\tsourceID\n";  
            print $O1 "SNP\tMapped_gene\tTrait_chr\tTrait_start\tTrait_end\tP\tSourceid\n";
            my $output  = join ("\n",@vs);
            print $O1 "$output\n";
            close $O1;
            my $pop = lc($Population);
            if ($pop =~/mix|eur/){ #混合人群用eur做
                $Population= "EUR";
                $pop = "eur";
                my $command= "plink --bfile /share/data0/1kg_phase3_v5_hg19/${Population}/1kg.phase3.v5.shapeit2.${pop}.hg19.all.SNPs.uniq_posID  --threads 20 --clump-r2 0.5 --clump-kb 250 --clump $fo1 --out ../output/tissue_pop_clump_0.5/$used_to_clump";
                # print "$command\n";    
                system $command;
                print "$fo1\n";                
            }
            else{
                my $command= "plink --bfile /share/data0/1kg_phase3_v5_hg19/${Population}/1kg.phase3.v5.shapeit2.${pop}.hg19.all.SNPs.uniq.posID  --threads 20 --clump-r2 0.5 --clump-kb 250 --clump $fo1 --out ../output/tissue_pop_clump_0.5/$used_to_clump";
                # print "$command\n";
                system $command;
                print "$fo1\n";
            }           
        }
    }
}