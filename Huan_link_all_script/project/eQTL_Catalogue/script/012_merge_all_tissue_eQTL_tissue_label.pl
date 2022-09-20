#合并部分数据集位点并利用"/share/Projects/huanhuan/ref_data/1kg_phase3_hg38/EUR/ID_MAF/snp/EUR_id_maf_more_than0.05.txt.gz"筛选落在eur maf >0.05的eqtl,同时利用"/share/Projects/huanhuan/ref_data/1kg_phase3_hg38/EUR/ID_MAF/snp/EUR_id_maf_greater0.txt.gz" 补全得"../output/01_merge_all_tissue_cis_eQTL_1kg_Completion.txt.gz"
#并得进入筛选的数据集信息"../output/need_study_for_hotspot_download_tabix_ftp_paths.tsv"
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;
use Parallel::ForkManager;

my $f1 = "../output/need_study_for_hotspot_download_tabix_ftp_paths.tsv";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 

my $fo2 = "../output/01_merge_all_tissue_cis_eQTLtissue_label.bed.gz";
open my $O2, "| gzip >$fo2" or die $!;
# print $O2 "SNP_chr\tstart\tend\tpvalue\ttissue\n";

my $dir = "/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/eQTL_Catalogue/data";
my (%hash1,%hash2,%hash3,%hash4);




my $f3 = "/share/Projects/huanhuan/ref_data/1kg_phase3_hg38/EUR/ID_MAF/snp/EUR_id_maf_more_than0.05.txt.gz";  #为了控制EUR人群
# open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
open( my $I3 ,"gzip -dc $f3|") or die ("can not open input file '$f3' \n"); #读压缩文件

while(<$I3>)
{
    chomp;
    unless(/^chr/){
        my @f = split/\t/;
        my $id =$f[0];
        my $maf =$f[1]; 
        my @t=split/:/,$id;
        my $variant =join("\t",@t);
        my $position = join("\t",@t[0..1]);
        $hash3{$variant}=1;
        $hash4{$position}=1;
    }
}



while(<$I1>)
{
    chomp;
    unless(/^study/){
        my @f = split/\t/;
        my $study= $f[0];
        my $qtl_group =$f[1];
        my $tissue_label = $f[4];
        my $condition_label =$f[5];
        my $quant_method =$f[6];
        my $ftp_path =$f[-2];
        my @t= split/\//,$ftp_path;
        my $filename=$t[-1];
        print "$study\t$filename\n";
        my $f2 = "${dir}/${study}/${filename}";
        open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件
        while(<$I2>)
        {
            chomp;
            unless(/^molecular_trait_id/){
                my @ss = split/\t/;
                my $chr = $ss[1];
                my $pos= $ss[2];
                my $ref= $ss[3];
                my $alt= $ss[4];
                # my $variant =$ss[5];
                my $pvalue =$ss[8];
                my $chrpos ="$chr\t$pos";
                my $variant = "$chrpos\t$ref\t$alt";
                my $start = $pos-1+1;
                my $end =$pos+1;
                my $bed="chr${chr}\t$start\t$end\t$pvalue\t$qtl_group";
                if($pvalue <5e-8){            
                    # print "$variant\n$chrpos\n";
                    if(exists $hash3{$variant}){  #test in eur 
                        print $O2 "$bed\n";

                    }
                    else{
                        if(exists $hash4{$chrpos}){
                            print $O2 "$bed\n";                            
                        }
                    }
                }
            }
        }
    }
}


close($O2);
print "start sort\n";
system "zless $fo2 |sort -k1,1 -k2,2n |gzip >../output/01_merge_all_tissue_cis_eQTLtissue_label_sorted.bed.gz";