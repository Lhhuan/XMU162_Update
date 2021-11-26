#!/usr/bin/perl
use strict;
use warnings;
my (@arr,$len,$a,$sample1,$sample2);
my $i=0;
#my $input=$ARGV[0];
my $out1=$ARGV[0]; # 外面传参的第一个参数
my $out2=$ARGV[1];# 外面传参的第二个参数
open (FIN,"gunzip -c /public/home/fengyang/data11/data11_all_cfy.cns.dp.gt.gz|") or die $!;
#open (FIN,"gunzip -c /public/home/fengyang/c.gz|") or die $!;

open (OUT,"|gzip >/public/home/fengyang/data11/$out1\.gz") or die $!;
open (OUT1,"|gzip >/public/home/fengyang/data11/$out2\.gz") or die $!;
while(<FIN>){
        chomp;
        @arr =split /\t/,$_;
        #$len =$#arr;
        $len =scalar(@arr);
        #print "$len\n";
        print OUT "$arr[0]\t";
        print OUT1 "$arr[0]\t";
        #print "$len\n";
        for $i(1..11278){
        $a=$arr[$i];
        if($a =~ /(\S+)\/(\S+)/){
                $sample1=$1;
                $sample2=$2;
                print OUT "$sample1\t";
                print OUT1 "$sample2\t";
}else{
                print OUT "\t";
                print OUT1 "\t";
}


}
        print OUT "\n";
        print OUT1 "\n";
}
close(FIN);
close(OUT1);
close(OUT2);