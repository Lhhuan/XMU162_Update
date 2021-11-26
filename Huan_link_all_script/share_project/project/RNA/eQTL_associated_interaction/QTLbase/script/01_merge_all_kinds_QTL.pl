#将../data/中的各种QTL merge到一起，得../output/01_all_kinds_QTL.txt
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;

my $fo1 = "../output/01_all_kinds_QTL1.txt.gz";
# open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
open my $O1, "| gzip >$fo1" or die $!;
print  $O1 "SNP_chr\tSNP_pos\tMapped_gene\tTrait_chr\tTrait_start\tTrait_end\tPvalue\tSourceid\tQTL_type\n";

#-------------------------------------------------------------#获取组织名称
my $dir = "/share/data0/QTLbase/data/";
opendir (DIR, $dir) or die "can't open the directory!";
my @files = readdir DIR; #获取一个文件夹下的所有文件
my @tissues;
my $suffix = "QTL.txt.gz"; #文件名后缀QTL.txt
foreach my $file(@files){
    if ($file =~/$suffix/){ 
        my $f1 = "/share/data0/QTLbase/data/$file";
        # open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n";
        open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
        while(<$I1>)
        {
            chomp;
            unless(/SNP_chr/){
                # print "$file\n";
                my @f =split/\t/;
                my $type = $file; 
                $type =~ s/\.txt\.gz//g;
                print $O1 "$_\t$type\n";
            }
        }
    }
}