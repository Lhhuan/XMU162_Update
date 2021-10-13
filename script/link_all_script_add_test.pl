#!/usr/bin/perl
use warnings;
use strict;
use utf8;
use File::Basename;

my $dir_out   = "/home/huanhuan/Script_backup/Huan_link_all_script";
mkdir $dir_out unless -d $dir_out;

chdir "/home/huanhuan/";
system "find -name *.pl > /home/huanhuan/perl_script";
system "find -name *.R > /home/huanhuan/R_script";
system "find -name *.sh > /home/huanhuan/all_sh_script";
system "find -name *.py > /home/huanhuan/python_script";
system "find -name *readme.txt > /home/huanhuan/readme";
system "find -name *.job > /home/huanhuan/job_script";

my $f1 ="/home/huanhuan/perl_script";
open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n";
while(<$I1>)
{
    chomp;
    my $script = $_;
    unless($script=~/^\.\/anaconda/){
        unless($script=~/^\.\/Script_backup/){
            unless($script=~/^\.\/php/){
                unless($script=~/^\.\/tools/){
                    unless($script=~ /^\.\/\.local/){
                        unless($script=~ /^\.\/R/){
                            my $file = basename($script);
                            my $dir = dirname($script);
                            $dir=~s/^\.\///;
                            my $do = "$dir_out/$dir";
                            # mkdir  $do unless -d $do;
                            unless(-e $do ){
                                system "mkdir -p $do";
                            }
                            my $new_file ="$do/$file";
                            unless(-e "$new_file"){
                                #print "#$new_file\n";
                                my $link1 = "ln \"$script\" \"$new_file\"" ;  #把变量引起来，这样就可以将名为12 3.pl的脚本copy 过来。而不会只copy12 而不copy 12 3.pl
                            system "$link1\n";
                                # print "$link1\n";
                            }
                        }
                    }
                }
            }
        }
    }
}

my $f2 ="/home/huanhuan/R_script";
open my $I2, '<', $f2 or die "$0 : failed to open input file '$f2' : $!\n";
while(<$I2>)
{
   chomp;
    my $script = $_;
    unless($script=~/^\.\/anaconda/){
        unless($script=~/^\.\/Script_backup/){
            unless($script=~/^\.\/php/){
                unless($script=~/^\.\/tools/){
                    unless($script=~ /^\.\/\.local/){
                        unless($script=~ /^\.\/R/){
                            my $file = basename($script);
                            my $dir = dirname($script);
                            $dir=~s/^\.\///;
                            my $do = "$dir_out/$dir";
                            # mkdir  $do unless -d $do;
                            unless(-e $do ){
                                system "mkdir -p $do";
                            }
                            my $new_file ="$do/$file";
                            unless(-e "$new_file"){
                                #print "#$new_file\n";
                                my $link1 = "ln \"$script\" \"$new_file\"" ;  #把变量引起来，这样就可以将名为12 3.pl的脚本copy 过来。而不会只copy12 而不copy 12 3.pl
                            system "$link1\n";
                                # print "$link1\n";
                            }
                        }
                    }
                }
            }
        }
    }
}

my $f3 ="/home/huanhuan/all_sh_script";
open my $I3, '<', $f3 or die "$0 : failed to open input file '$f3' : $!\n";
while(<$I3>)
{
    chomp;
    my $script = $_;
    unless($script=~/^\.\/anaconda/){
        unless($script=~/^\.\/Script_backup/){
            unless($script=~/^\.\/php/){
                unless($script=~/^\.\/tools/){
                    unless($script=~ /^\.\/\.local/){
                        unless($script=~ /^\.\/R/){
                            my $file = basename($script);
                            my $dir = dirname($script);
                            $dir=~s/^\.\///;
                            my $do = "$dir_out/$dir";
                            # mkdir  $do unless -d $do;
                            unless(-e $do ){
                                system "mkdir -p $do";
                            }
                            my $new_file ="$do/$file";
                            unless(-e "$new_file"){
                                #print "#$new_file\n";
                                my $link1 = "ln \"$script\" \"$new_file\"" ;  #把变量引起来，这样就可以将名为12 3.pl的脚本copy 过来。而不会只copy12 而不copy 12 3.pl
                            system "$link1\n";
                                # print "$link1\n";
                            }
                        }
                    }
                }
            }
        }
    }
}

my $f4 ="/home/huanhuan/readme";
open my $I4, '<', $f4 or die "$0 : failed to open input file '$f4' : $!\n";
while(<$I4>)
{
    chomp;
    my $script = $_;
    unless($script=~/^\.\/anaconda/){
        unless($script=~/^\.\/Script_backup/){
            unless($script=~/^\.\/php/){
                unless($script=~/^\.\/tools/){
                    unless($script=~ /^\.\/\.local/){
                        unless($script=~ /^\.\/R/){
                            my $file = basename($script);
                            my $dir = dirname($script);
                            $dir=~s/^\.\///;
                            my $do = "$dir_out/$dir";
                            # mkdir  $do unless -d $do;
                            unless(-e $do ){
                                system "mkdir -p $do";
                            }
                            my $new_file ="$do/$file";
                            unless(-e "$new_file"){
                                #print "#$new_file\n";
                                my $link1 = "ln \"$script\" \"$new_file\"" ;  #把变量引起来，这样就可以将名为12 3.pl的脚本copy 过来。而不会只copy12 而不copy 12 3.pl
                            system "$link1\n";
                                # print "$link1\n";
                            }
                        }
                    }
                }
            }
        }
    }
}

my $f5 ="/home/huanhuan/python_script";
open my $I5, '<', $f5 or die "$0 : failed to open input file '$f5' : $!\n";
while(<$I5>)
{
    chomp;
    my $script = $_;
    unless($script=~/^\.\/anaconda/){
        unless($script=~/^\.\/Script_backup/){
            unless($script=~/^\.\/php/){
                unless($script=~/^\.\/tools/){
                    unless($script=~ /^\.\/\.local/){
                        unless($script=~ /^\.\/R/){
                            my $file = basename($script);
                            my $dir = dirname($script);
                            $dir=~s/^\.\///;
                            my $do = "$dir_out/$dir";
                            # mkdir  $do unless -d $do;
                            unless(-e $do ){
                                system "mkdir -p $do";
                            }
                            my $new_file ="$do/$file";
                            unless(-e "$new_file"){
                                #print "#$new_file\n";
                                my $link1 = "ln \"$script\" \"$new_file\"" ;  #把变量引起来，这样就可以将名为12 3.pl的脚本copy 过来。而不会只copy12 而不copy 12 3.pl
                            system "$link1\n";
                                # print "$link1\n";
                            }
                        }
                    }
                }
            }
        }
    }
}



my $f6 ="/home/huanhuan/job_script";
open my $I6, '<', $f6 or die "$0 : failed to open input file '$f6' : $!\n";
while(<$I6>)
{
    chomp;
    my $script = $_;
    unless($script=~/^\.\/anaconda/){
        unless($script=~/^\.\/Script_backup/){
            unless($script=~/^\.\/php/){
                unless($script=~/^\.\/tools/){
                    unless($script=~ /^\.\/\.local/){
                        unless($script=~ /^\.\/R/){
                            my $file = basename($script);
                            my $dir = dirname($script);
                            $dir=~s/^\.\///;
                            my $do = "$dir_out/$dir";
                            # mkdir  $do unless -d $do;
                            unless(-e $do ){
                                system "mkdir -p $do";
                            }
                            my $new_file ="$do/$file";
                            unless(-e "$new_file"){
                                #print "#$new_file\n";
                                my $link1 = "ln \"$script\" \"$new_file\"" ;  #把变量引起来，这样就可以将名为12 3.pl的脚本copy 过来。而不会只copy12 而不copy 12 3.pl
                            system "$link1\n";
                                # print "$link1\n";
                            }
                        }
                    }
                }
            }
        }
    }
}