#!/usr/bin/perl
use strict;
use warnings;

binmode STDOUT, ":utf8";
use utf8;

use JSON;

# my $json ={
#     'book' => {
#         'title' => 'smth',
#         'num_page' => 234
#     }
# };

# my $data = decode_json($json);
# open my $fh, ">", "data_out.json";
# print $fh encode_json($data);
# close $fh;


# my $json ={
#     'book' => [
#         {'title' => 'smth'},
#         {'num_page' => 234}
#     ]
# };
# open my $fh, ">", "data_out.json";
# print $fh encode_json($json);
# close $fh;

my $name = "test";
my $type = "A";
my $data = "1.1.1.1";
my $ttl  = 84600;

my %rec_hash = ('name'=>$name, 'type'=>$type,'data'=>$data,{'ttl'=>$ttl,'ttl'=>"500"});
my %hash1;
push @{$hash1{'name'}},$name;
push @{$hash1{'type'}},$type;
push @{$hash1{'data'}},$data;
push @{$hash1{'ttl'}},84600;
push @{$hash1{'ttl'}},$ttl;
# my $json = encode_json \%rec_hash;
# print $json;
# # print "%rec_hash\n";


foreach my $k(sort keys %hash1){
    my @vs =@{$hash1{$k}};
    my $v =join("\t",@vs);
    print "$k\t$v\n";
}
my $json = encode_json \%hash1;
print $json;
