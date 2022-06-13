#通过"/home/huanhuan/project/link_database/JEME/ENCODE_Roadmap_info_man.txt" 筛选normal tissue/cell 将"/home/huanhuan/project/link_database/OncoBase/EpiTensor.txt.gz" 中的TSS_TSS等type分开，得"./output/01_${type}.bed.gz，转hg38得"./output/hg38/01_${type}.bed",sort得"./output/hg38/01_${type}_sorted.bed.gz"
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;
use Parallel::ForkManager;

my (%hash1,%hash2);

my $f2 = "/home/huanhuan/project/link_database/JEME/ENCODE_Roadmap_info_man.txt";
open my $I2, '<', $f2 or die "$0 : failed to open input file '$f2' : $!\n"; 

while(<$I2>)
{
    chomp;
    unless(/^ID/){
        my @f = split/\t/;
        for (my $i=0;$i<8;$i++){ #对文件进行处理，把所有未定义的空格等都替换成NONE
            unless(defined $f[$i]){$f[$i] = "NA"};
        }
        my $ID = $f[0];
        my $Sample_name =$f[1];
        my $cancer =$f[-1];
        unless($cancer =~/YES/){
            $hash1{$Sample_name}=1;
            # print "$ID\n";
        }
    }
}

my @types=("TSS_TSS","ENH_ENH","TSS_ENH");
foreach my $type(@types){
    my $f1 = "/home/huanhuan/project/link_database/OncoBase/EpiTensor.txt.gz";
    open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
    # open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
    my $fo1 = "./output/01_${type}.bed.gz";
    # open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";
    open my $O1, "| gzip >$fo1" or die $!;
    while(<$I1>)
    {
        chomp;
        my @f=split/\t/;
        my $pos1= join("\t",@f[0..2]);
        my $pos2=join("\t",@f[3..5]);
        my $type1= $f[6];
        my $tissue =$f[7];
        my @t=split/_/,$tissue;
        my $ID= $t[0];
        
        if(exists $hash1{$ID}){
            # print "$ID\n";
            my $gene1=$f[8];
            my $gene2 =$f[9];
            if($type1 eq $type){
            # if($type1 =~/TSS_ENH/){
                print $O1 "$pos1\t$gene2\t$tissue\n";
                print $O1 "$pos2\t$gene1\t$tissue\n";
            }
        }
    }
    close($O1);
    my $f2="./output/hg38/01_${type}.bed";
    my $f3 = "./output/hg38/01_${type}_sorted.bed.gz";
    # my $f4= ""
    system "liftOver $fo1  /home/huanhuan/reference/hg19ToHg38.over.chain.gz $f2  ./output/01_unmap_${type}_enhancer_target_sample.bed" ;
    system "zless $f2 |sort  -k1,1 -k2,2n|gzip >$f3";
    print "$type finish\n";
}
