#download getx related marks to "/share/data0/GTEx/annotation/ROADMAP/sample/${roadmap_ID}
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use List::Util qw/sum/;
use Env qw(PATH);
use Parallel::ForkManager;

my $f1 = "../tissue_roadmapID.txt";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
# open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件

my $fo1 = "unique_roadmap_id.txt";
open my $O1, '>', $fo1 or die "$0 : failed to open output file '$fo1' : $!\n";

my (%hash1,%hash2,%hash3,%hash4);


my @markers = ("H3K4me1","H3K4me3","H3K9ac","H3K9me3","H3K27ac","H3K27me3","H3K36me3");
# my @markers = ("H3K4me1");

my $pm = Parallel::ForkManager->new(10); ## 设置最大的线程数目
while(<$I1>)
{
    chomp;
    unless(/^GTEx/){
        my @f=split/\t/;
        my $getx_tissue = $f[0];
        my $roadmap_tissue_ID = $f[1];
        my @ts =  split/,/,$roadmap_tissue_ID;
        unless($roadmap_tissue_ID   =~/\bNA\b/){
            foreach my $t(@ts){
                $t =~ s/\(/;/g;
                $t =~ s/\)/;/g;
                my @tts = split/;/,$t;
                my $roadmap_ID =$tts[1];
                unless(exists $hash1{$roadmap_ID}){
                    $hash1{$roadmap_ID}=1;
                    print $O1 "$roadmap_ID\n";
                    my $output_dir = "/share/data0/GTEx/annotation/ROADMAP/sample/${roadmap_ID}";
                    if(-e $output_dir){
                        print "${output_dir}\texist\n";
                    }
                    else{
                        system "mkdir -p $output_dir";
                    }
                    $ENV{'output_dir'} = $output_dir ;
                    $ENV{'roadmap_ID'} = $roadmap_ID ;
                    print "$roadmap_ID\n";
                    chdir $output_dir;
                    
                    foreach my $marker(@markers){
                        my $pid = $pm->start and next; #开始多线程
                        $ENV{'marker'} = $marker ;
                        system "bash /home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/ROADMAP/script/download_markers.sh";
                        print "$marker\n";
                        $pm->finish;  #多线程结束
                    }
                }
            }
        }
    }
}