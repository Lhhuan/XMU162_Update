#../../output/five_race_1kg_snp.txt.gz  补全"${dir}/${tissue}${suffix}"; 得"../../output/${tissue}_trans_eQTL_1kg_Completion.txt.gz"
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;
use Parallel::ForkManager;

# my @tissues =  ("ACC","BRCA","COAD","ESCA","KICH","KIRC","KIRP","LAML","LIHC","LUAD","LUSC","OV","PAAD","PRAD","SKCM","STAD","TGCT","THCA","UCEC","UCS");
my @tissues =  ("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM");



my $pm = Parallel::ForkManager->new(5);
foreach my $tissue(@tissues){
    my $pid = $pm->start and next; #开始多线程
    # my $tissue = "Whole_Blood";
    print "$tissue\tstart\n";
    my $f1 ="/home/huanhuan/project/RNA/eQTL_associated_interaction/PancanQTL/data/trans_eQTLs_all_re.gz";
    open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件

    my $f2 = "../../output/five_race_1kg_snp.txt.gz";
    # open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
    open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件

    # $tissue =~s/-/_/g;
    my $output_dir = "../../output/${tissue}";
    unless(-e $output_dir){
        system "mkdir -p $output_dir";
    }
    #---------------
    my $fo1 = "$output_dir/${tissue}_trans_eQTL_1kg_Completion.txt.gz";
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
        unless(/^SNP_chr/){
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

