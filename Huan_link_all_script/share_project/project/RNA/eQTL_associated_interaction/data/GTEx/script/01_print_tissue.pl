#将../raw_data/GTEx_Analysis_v8_eQTL/*.v8.signif_variant_gene_pairs.txt.gz进行初步normal, 得gene name 和部分RSid从对应的v8.egenes.txt.gz中获得，得文件../output/01_normal_format.txt.gz
#并得unique的varint POS的vcf文件../output/01_unique_variant-id_position_hg38.vcf.gz
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;


#-------------------------------------------------------------#获取组织名称
my $dir = "../raw_data/GTEx_Analysis_v8_eQTL";
opendir (DIR, $dir) or die "can't open the directory!";
my @files = readdir DIR; #获取一个文件夹下的所有文件
my @tissues;
my $suffix = ".v8.signif_variant_gene_pairs.txt.gz";
foreach my $file(@files){
    my $suffix = ".v8.signif_variant_gene_pairs.txt.gz";
    if ($file =~/$suffix/){ 
        $file =~ s/$suffix//g; #提取前缀disease
        my $tissue =$file;
        print "$tissue\n";
        push @tissues,$tissue;
    }
}