#!/usr/bin/perl
# use strict;
use warnings;
my @inter = (6,7,12,13);
our (%hash,%hash6,%hash7,%hash12,%hash13);

foreach my $j(@inter){
	
	$hash{$j}{"kkk"}=1;
	# print "$hash{$j}{kkk}\n";
	print $hash6{"kkk"};
	# print "$j\n";

}

#}

