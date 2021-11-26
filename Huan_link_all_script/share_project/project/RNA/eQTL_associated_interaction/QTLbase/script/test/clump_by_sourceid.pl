#分source进行clump,将不同类型QTL按照sourceid分割成"../output/used_to_clump/$used_to_clump"，然后用对应的人种文件/share/data0/1kg_phase3_v5_hg19/${Population}/1kg.phase3.v5.shapeit2.${pop}.hg19.all.SNPs.uniq.posID
#进行clump, 得../output/clump/${used_to_clump}.clump

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
        my $Population =$f[4];
        $hash1{$Sourceid}=$Population;
    }
}

#-------------------------------------------------------------#获取组织名称
my $dir = "/share/data0/QTLbase/data/";
opendir (DIR, $dir) or die "can't open the directory!";
my @files = readdir DIR; #获取一个文件夹下的所有文件
my @tissues;
my $suffix = "QTL.txt.gz"; #文件名后缀QTL.txt
foreach my $file(@files){
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
                my $Pvalue =$f[-2];
                my $Sourceid =$f[-1];
                my $SNP = "${SNP_chr}_${SNP_pos}";
                my $P = $Pvalue;
                my $v = "$SNP\t$P";
                push @{$hash2{$Sourceid}},$v;
            }
        }
        my $type = $file; 
        $type =~ s/\.txt.gz//g;
        foreach my $id(sort keys %hash2){
            #------------------------------------------
            my $used_to_clump = "${type}_${id}";
            my $fo1 = "../output/used_to_clump/$used_to_clump";
            open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
            print $O1 "SNP\tP\n";
            #------------------------------------------------
            my @vs = @{$hash2{$id}};
            my $output  = join ("\n",@vs);
            print $O1 "$output\n";
            close $O1;
            if (exists $hash1{$id}){ #用对应人群的1KG数据进行clump

                my $Population = $hash1{$id};
                my $pop = lc($Population);
                if ($pop =~/mix|eur/){ #混合人群用eur做
                    $Population= "EUR";
                    $pop = "eur";
                    my $command= "plink --bfile /share/data0/1kg_phase3_v5_hg19/${Population}/1kg.phase3.v5.shapeit2.${pop}.hg19.all.SNPs.uniq_posID  --clump-r2 0.2 --clump-kb 250 --clump $fo1 --out ../output/clump/$used_to_clump";
                    # print "$command\n";    
                    system $command;
                    print "$fo1\n";                
                }
                else{
                    my $command= "plink --bfile /share/data0/1kg_phase3_v5_hg19/${Population}/1kg.phase3.v5.shapeit2.${pop}.hg19.all.SNPs.uniq.posID  --clump-r2 0.2 --clump-kb 250 --clump $fo1 --out ../output/clump/$used_to_clump";
                    # print "$command\n";
                    system $command;
                    print "$fo1\n";
                }
            }
        }
    }
}