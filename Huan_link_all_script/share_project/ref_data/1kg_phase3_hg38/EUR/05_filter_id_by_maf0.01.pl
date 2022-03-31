#将"./ID_MAF/EUR_2020-08-05_chr${i}.filtered.shapeit2-duohmm-phased.maf.id.txt.gz"  中的EUR提出来，得EUR_sample_list.txt

#!/usr/bin/perl
use warnings;
use strict; 
use utf8;

my $fo1 = "./ID_MAF/snp/EUR_id_maf_more_than0.01.txt.gz";
open my $O1, "| gzip >$fo1" or die $!;
my (%hash1,%hash2);

print $O1 "chr:pos:ref:alt\tmaf\trsid\n";


my @numbers=(1..22);
foreach my $i(@numbers){
    print "$i\n";
    my $f2 = "./ID_MAF/EUR_2020-08-05_chr${i}.filtered.shapeit2-duohmm-phased.maf.id.txt.gz";
    open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件
    while(<$I2>)
    {
        chomp;
        unless(/^#/){
            my @f =split/\t/;
            my $id =$f[0];
            my $maf =$f[1];
            # print "$_\n";
            if($maf >=0.01){ 
                print $O1 "$_\n";
                # if(exists $hash1{$id}){
                #     my $rsid =$hash1{$id};
                #     print $O1 "$_\t$rsid\n";
                # }
            }
        }
    }
}

