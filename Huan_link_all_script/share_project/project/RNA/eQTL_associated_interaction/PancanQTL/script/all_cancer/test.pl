#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use File::Basename;
use List::Util qw/max min/;
use Env qw(PATH);
use Parallel::ForkManager;
use List::MoreUtils ':all';


my %hash1;
$hash1{"aa"}=1;
$hash1{"nb"}=1;
my $hash_length = keys %hash1;
print "$hash_length\n";


# my @list=qw /1 2 3 2 1 4 aa a bb c  b bb d/;
# foreach (@list){print "$_ ";}

# print "\n###################\n";
# my @uni=uniq(@list);
# foreach (@uni){print "$_ ";}
my $f2 = "${out_dir}/${marker}_${input_file_base_name}";
my %hash2;
                while(<$I2>)
                {
                    chomp;
                    my @f = split/\t/;
                    my $hotspot_chr = $f[0];
                    my $hotspot_start = $f[1];
                    my $hotspot_end = $f[2];
                    my $factor_chr = $f[3];
                    my $factor_start = $f[4];
                    my $factor_end = $f[5];
                    my $overlap_bp = $f[-1];
                    my $factor_length = $factor_end - $factor_start;
                    my $k = "$hotspot_chr\t$hotspot_start\t$hotspot_end";
                    my $v = "$factor_chr\t$factor_start\t$factor_end\t$overlap_bp";
                    push @{$hash2{$k}},$v;
                }
                foreach my $k(sort keys %hash2){
                    print "$k\n";
                }