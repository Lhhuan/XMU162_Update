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


my @files = ("/home/huanhuan/perl_script","/home/huanhuan/R_script","/home/huanhuan/all_sh_script","/home/huanhuan/python_script","/home/huanhuan/readme","/home/huanhuan/job_script");
foreach my $f1(@files){
# my $f1 ="/home/huanhuan/perl_script";
    open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n";
    while(<$I1>)
    {
        chomp;
        my $script = $_;
        unless($script=~/^\.\/anaconda|^\.\/Script_backup|^\.\/php|^\.\/tools|^\.\/\.local|^\.\/R|^\.\/\.conda|^\.\/\.config|^\.\/\.aspera/){  
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

#-------------------------------------------share before
my $dir_out2   = "/home/huanhuan/Script_backup/Huan_link_all_script/share_project";
mkdir $dir_out2 unless -d $dir_out2;

chdir "/share/Projects/huanhuan/";
system "find -name *.pl > /home/huanhuan/share_script";
system "find -name *.R >> /home/huanhuan/share_script";
system "find -name *.sh >> /home/huanhuan/share_script";
system "find -name *.py >> /home/huanhuan/share_script";
system "find -name *readme.txt >> /home/huanhuan/share_script";
system "find -name *.job >> /home/huanhuan/share_script";

my $f1 ="/home/huanhuan/share_script";
    open my $I1, '<', $f1 or die "$0 : failed to open input file '$f1' : $!\n";
    while(<$I1>)
    {
        chomp;
        my $script = $_;
        unless($script=~/^\.\/anaconda|^\.\/Script_backup|^\.\/php|^\.\/tools|^\.\/\.local|^\.\/R|^\.\/\.conda|^\.\/miniconda3/){  
            my $file = basename($script);
            my $dir = dirname($script);
            $dir=~s/^\.\///;
            my $do = "$dir_out2/$dir";
            # mkdir  $do unless -d $do;
            unless(-e $do ){
                system "mkdir -p $do";
            }
            my $new_file ="$do/$file";
            unless(-e "$new_file"){
                #print "#$new_file\n";
                my $link1 = "cp \"$script\" \"$new_file\"" ;  #把变量引起来，这样就可以将名为12 3.pl的脚本copy 过来。而不会只copy12 而不copy 12 3.pl
            system "$link1\n";
                # print "$link1\n";
            }
        }
    }

