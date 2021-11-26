#将"../../../output/ALL_eQTL/cis_trans/QTLbase_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_cis_trans_eQTL.txt.gz" 处理成cis_trans_1MB, cis_trans_10MB的格式，
#得../../../output/ALL_eQTL/cis_trans/QTLbase_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_cis_trans_eQTL_normalized.txt.gz
#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use List::Util qw/max min/;
use 5.010;


my $f1 = "../../../output/ALL_eQTL/cis_trans/QTLbase_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_cis_trans_eQTL.txt.gz";
# open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n"; 
open( my $I1 ,"gzip -dc $f1|") or die ("can not open input file '$f1' \n"); #读压缩文件
my $fo1 = "../../../output/ALL_eQTL/cis_trans/QTLbase_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_cis_trans_eQTL_normalized.txt.gz";
open my $O1, "| gzip >$fo1" or die $!;
my $fo2 = "../../../output/ALL_eQTL/cis_trans/value_more_than1.txt";
# open my $O2, "| gzip >$fo2" or die $!;
open my $O2, '>', $fo2 or die "$0 : failed to open output file '$fo2' : $!\n";
my %hash1;
print $O1 "chr\tpos\t1MB_emplambda\tcis_trans_1MB\t10MB_emplambda\tcis_trans_10MB\n";
# my %name;
while(<$I1>)
{
    chomp;
    unless(/^emplambda/){
        my @f = split/\t/;
        my $emplambda =$f[0];
        my $pos =$f[1];
        my $chr =$f[2];
        my $cis_or_trans = $f[3];
        my $k = "$chr\t$pos";
        my $v = "$emplambda\t$cis_or_trans";
        push @{$hash1{$k}},$v;
    }
}


foreach my $k (sort keys %hash1){
    my @vs = @{$hash1{$k}};
    my %hash2;
    @vs = grep { ++$hash2{$_} < 2 } @vs;
    print $O1 "$k\t";
    my @m1s= ();
    my @m10s =();
    foreach my $v(@vs){
        my @t = split/\t/,$v;
        my $emplambda = $t[0];
        if($v =~ /cis_1MB/){
            my $M1 = "$emplambda\tcis";
            push @m1s,$M1;
        }
        elsif($v =~ /trans_1MB/){
            my $M1 = "$emplambda\ttrans";
            push @m1s,$M1;
        }
        elsif($v =~ /cis_10MB/){
            my $M10 = "$emplambda\tcis";
            push @m10s,$M10;
        }
        else{
            my $M10 = "$emplambda\ttrans";
            push @m10s,$M10;            
        }
    }
    my $m1_number =@m1s;
    my $m10_number =@m10s;
    if ($m1_number <2 && $m10_number <2){
        my $m1 = $m1s[0];
        my $m2 = $m10s[0];
        print $O1 "$k\t$m1\t$m2\n";
    }
    else{
        print $O2 "$k\t$m1_number\t$m10_number\n";
        print $O2 "$k\t@m1s\t@m10s\n";
    }

}