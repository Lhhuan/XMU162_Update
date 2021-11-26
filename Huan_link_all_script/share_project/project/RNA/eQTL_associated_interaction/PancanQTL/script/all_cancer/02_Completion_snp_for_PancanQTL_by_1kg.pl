 # 用"/share/data0/1kg_phase3_v5_hg19/EUR/1kg.phase3.v5.shapeit2.eur.hg19.all.SNPs.vcf.gz" 补全"${dir}/${tissue}${suffix}"; 得"../../output/${tissue}_cis_eQTL_1kg_Completion.txt.gz"
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;
use Parallel::ForkManager;

my @tissues =  ("ACC","BRCA","COAD","ESCA","KICH","KIRC","KIRP","LAML","LIHC","LUAD","LUSC","OV","PAAD","PRAD","SKCM","STAD","TGCT","THCA","UCEC","UCS");

my $pm = Parallel::ForkManager->new(5);
foreach my $tissue(@tissues){
    my $pid = $pm->start and next; #开始多线程
    # my $tissue = "Whole_Blood";
    print "$tissue\tstart\n";
    my $f1 ="/home/huanhuan/project/RNA/eQTL_associated_interaction/PancanQTL/data/cis_eQTLs_all_re.gz";
    open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件

    my $f2 = "/share/data0/1kg_phase3_v5_hg19/EUR/1kg.phase3.v5.shapeit2.eur.hg19.all.SNPs.vcf.gz";
    # open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
    open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件

    my $f3 = "/share/data0/1kg_phase3_v5_hg19/EAS/1kg.phase3.v5.shapeit2.eas.hg19.all.SNPs.vcf.gz";
    open( my $I3 ,"gzip -dc $f3|") or die ("can not open input file '$f3' \n"); #读压缩文件

    my $f4 = "/share/data0/1kg_phase3_v5_hg19/SAS/1kg.phase3.v5.shapeit2.sas.hg19.all.SNPs.vcf.gz";
    open( my $I4 ,"gzip -dc $f4|") or die ("can not open input file '$f4' \n"); #读压缩文件

    my $f5 = "/share/data0/1kg_phase3_v5_hg19/AMR/1kg.phase3.v5.shapeit2.amr.hg19.all.SNPs.vcf.gz";
    open( my $I5 ,"gzip -dc $f5|") or die ("can not open input file '$f5' \n"); #读压缩文件

    my $f6 = "/share/data0/1kg_phase3_v5_hg19/AFR/1kg.phase3.v5.shapeit2.afr.hg19.all.SNPs.vcf.gz";
    open( my $I6 ,"gzip -dc $f6|") or die ("can not open input file '$f6' \n"); #读压缩文件
    #---------------------------mkdir
    # $tissue =~s/-/_/g;
    my $output_dir = "../../output/${tissue}";
    unless(-e $output_dir){
        system "mkdir -p $output_dir";
    }
    #---------------
    my $fo1 = "$output_dir/${tissue}_cis_eQTL_1kg_Completion.txt.gz";
    # open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
    open my $O1, "| gzip >$fo1" or die $!;
    print $O1 "SNP_chr\tSNP_pos\tPvalue\n";


    my %hash1;
    while(<$I1>)
    {
        chomp;
        unless(/^cancer_type/){
            my @f = split/\t/;
            my $cancer_type =$f[0];
            my $SNP_chr =$f[2];
            $SNP_chr =~ s/chr//g;
            my $SNP_pos =$f[3];
            my $Pvalue =$f[-1];
            if($cancer_type eq $tissue){
                # print "$cancer_type\n";
                my $k = "$SNP_chr\t$SNP_pos";
                $hash1{$k}=1;
                print $O1 "$k\t$Pvalue\n";
            # print "$k\t$Pvalue\n";
            }
        }
    }


    while(<$I2>)
    {
        chomp;
        unless(/^#/){
            my @f = split/\t/;
            my $CHROM =$f[0];
            my $POS =$f[1]; 
            my $pvalue = 0.05;
            my $k = "$CHROM\t$POS";
            unless (exists $hash1{$k}){
                print $O1 "$k\t$pvalue\n";
                $hash1{$k}=1;
            }
        }
    }

    while(<$I3>)
    {
        chomp;
        unless(/^#/){
            my @f = split/\t/;
            my $CHROM =$f[0];
            my $POS =$f[1]; 
            my $pvalue = 0.05;
            my $k = "$CHROM\t$POS";
            unless (exists $hash1{$k}){
                print $O1 "$k\t$pvalue\n";
                $hash1{$k}=1;
            }
        }
    }

    while(<$I4>)
    {
        chomp;
        unless(/^#/){
            my @f = split/\t/;
            my $CHROM =$f[0];
            my $POS =$f[1]; 
            my $pvalue = 0.05;
            my $k = "$CHROM\t$POS";
            unless (exists $hash1{$k}){
                print $O1 "$k\t$pvalue\n";
                $hash1{$k}=1;
            }
        }
    }

    while(<$I5>)
    {
        chomp;
        unless(/^#/){
            my @f = split/\t/;
            my $CHROM =$f[0];
            my $POS =$f[1]; 
            my $pvalue = 0.05;
            my $k = "$CHROM\t$POS";
            unless (exists $hash1{$k}){
                print $O1 "$k\t$pvalue\n";
                $hash1{$k}=1;
            }
        }
    }

    while(<$I6>)
    {
        chomp;
        unless(/^#/){
            my @f = split/\t/;
            my $CHROM =$f[0];
            my $POS =$f[1]; 
            my $pvalue = 0.05;
            my $k = "$CHROM\t$POS";
            unless (exists $hash1{$k}){
                print $O1 "$k\t$pvalue\n";
                # $hash1{$k}=1;
            }
        }
    }

    print "$tissue\tend\n";
    $pm->finish;  #多线程结束
}

