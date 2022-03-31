#利用"/share/data0/GTEx/data/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz" 将/share/data0/GTEx/data/GTEx_Analysis_v8_sQTL/*.v8.sqtl_signifpairs.txt.gz转换为hg19, 得/share/data0/GTEx/data/GTEx_Analysis_v8_sQTL_hg19/*.v8.sqtl_signifpairs.txt.gz
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use Parallel::ForkManager;

#-------------------------------------------------------------#获取组织名称
my $dir = "/share/data0/GTEx/data/GTEx_Analysis_v8_sQTL/";
opendir (DIR, $dir) or die "can't open the directory!";
my @files = readdir DIR; #获取一个文件夹下的所有文件
my @tissues;

foreach my $file(@files){
    my $suffix = ".v8.sqtl_signifpairs.txt.gz";
    if ($file =~/$suffix/){ 
        $file =~ s/$suffix//g; #提取前缀disease
        my $tissue =$file;
        push @tissues,$tissue;
    }
}
my $suffix = ".v8.sqtl_signifpairs.txt.gz";

# my $t =join("\n",@tissues);
# print "$t\n";
#---------------------------------------------------#对文件进行处理
my $f1 = "/share/data0/GTEx/data/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz";
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
my(%hash1,%hash2,%hash3,%hash4);


while(<$I1>){
    chomp;
    my @f = split/\t/;
    unless(/^variant_id/){
        my $variant_id_38 = $f[0];
        my $variant_id_37 = $f[-1];
        
        # print "$variant_id_38\t$variant_id_37\n";
        my $variant_id = $variant_id_37;
        $variant_id =~ s/chr//g;
        unless ($variant_id =~/NA/){
            # print "$variant_id\n";
            $hash1{$variant_id_38}=$variant_id_37;
            my @t = split/\_/,$variant_id_37;
            my $Chr =$t[0];
            my $Pos =$t[1];
            my $Ref=$t[2];
            my $Alt =$t[3];
            # print "$Chr\n";
        }
    }
}
my $pm = Parallel::ForkManager->new(20);

foreach my $tissue(@tissues){
    my $pid = $pm->start and next; #开始多线程
    # my $tissue = "Whole_Blood";
    my $f2 ="${dir}/${tissue}${suffix}";
    open( my $I2 ,"gzip -dc $f2|") or die ("can not open input file '$f2' \n"); #读压缩文件
    #-------------------
    my $fo1 = "/share/data0/GTEx/data/GTEx_Analysis_v8_sQTL_hg19/${tissue}${suffix}";
    open my $O1, "| gzip >$fo1" or die $!;
    my $fo2 = "/share/data0/GTEx/data/GTEx_Analysis_v8_sQTL_hg19/not_transfrom_from_38_to_19/${tissue}${suffix}";
    open my $O2, "| gzip >$fo2" or die $!;
    #---------------------
    while(<$I2>){
        chomp;
        my @f =split/\t/;
        my $Variant_id_38 =$f[0];
        my $info = join("\t",@f[1..$#f]);
        if (/^variant_id/){
            print $O1 "variant_id\tChr\tPos\t$info\n";
        }
        else{
            if (exists $hash1{$Variant_id_38}){
                my $variant_id_37 = $hash1{$Variant_id_38};
                # print "$variant_id_37\t$f2\n";
                my @t = split/\_/,$variant_id_37;
                my $Chr =$t[0];
                my $Pos =$t[1];
                my $Ref=$t[2];
                my $Alt =$t[3];
                print $O1 "$variant_id_37\t$Chr\t$Pos\t$info\n";
            }
            else{
                print $O2 "$Variant_id_38\n";
            }
        }
    }
    close $I2;
    close $O1;
    close $O2;
    $pm->finish;  #多线程结束
}

close $I1;