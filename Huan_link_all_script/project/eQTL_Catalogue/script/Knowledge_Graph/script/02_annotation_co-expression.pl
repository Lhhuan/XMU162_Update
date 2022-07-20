#为./output/01_hotspot_target_gene_reactomeFI.bed.gz  annotation co-expression数据，得./ENSG_G16808_S85825/${file_entrezgene}_${file_ensembl}，得./output/02_hotspot_target_gene_reactomeFI_co-expression.bed.gz 
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::MoreUtils ':all';

# my $f4 = "1234.bed.gz";
my $f4 = "../output/01_hotspot_target_gene_reactomeFI.bed.gz";
open( my $I4 ,"gzip -dc $f4|") or die ("can not open input file '$f4' \n"); #读压缩文件
my $fo1 = "../output/02_hotspot_target_gene_reactomeFI_co-expression.bed.gz"; #
open my $O1, "| gzip >$fo1" or die $!;
my %hash1;


my $dir = "/home/huanhuan/project/link_database/COXPRESdb/ENSG_G16808_S85825/";
opendir (DIR, $dir) or die "can't open the directory!";
my @dir = readdir DIR;
foreach my $file(@dir){
    if($file=~/same_pos/){
       my @t= split/_/,$file;
       my $ensg =$t[1];
       $ensg =~s/\.txt.*//g;
    #    print "$ensg\n";
       my $f1 = "$dir/$file";
    #    open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n";
       open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n");
       while(<$I1>)
        {
            chomp;
            my @f =split/\t/;
            my $co_ensg =$f[0];
            push @{$hash1{$ensg}},$co_ensg;
        }
    }
}

while(<$I4>)
{
    chomp;
    my @f= split/\t/;
    if(/^hotspot_chr/){
        print $O1 "$_\tegene_co_expression_gene\n";
        # print "$_\n";
    }
    else{
        my $egene = $f[3];
        for (my $i=4;$i<7;$i++){ #对文件进行处理，把所有未定义的空格等都替换成NONE
            unless(defined $f[$i]){$f[$i] = "NA"};
        }
        my $refine_line = join("\t",@f[0..6]);
        if(exists $hash1{$egene}){
            my @co_express =@{$hash1{$egene}};
            my $num= @co_express;
            if($num >0){
                my @co_exp =uniq(@co_express);
                # foreach my $aaa(@co_exp){
                #     print "$aaa\t";
                # }
                # print "\n";
                my $co_ensg = join(";",@co_exp);
                print $O1 "$refine_line\t$co_ensg\n";
            }
            else{
                print $O1 "$refine_line\tNA\n";
            }
        }
        else{
            print $O1 "$refine_line\tNA\n";
        }
    }
}