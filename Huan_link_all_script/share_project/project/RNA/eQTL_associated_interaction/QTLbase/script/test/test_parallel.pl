#!/usr/bin/perl
use warnings;
use strict; 
use utf8;
use Parallel::ForkManager; #多线程并行

use LWP::Simple;
# use threads;

# print localtime(time),"n";  #输出系统时间；
# my $j=0;
# my $thread;
# while()
# {    last if($j>=10)；这里控制一下任务数量，共10个；
#     #控制创建的线程数，这里是5，scalar函数返回列表threads->list()元素的个数；
#     while(scalar(threads->list())<5)  
#     {    $j++;
#          #创建一个线程，这个线程其实就是调用（引用）函数“ss”；
#          #函数‘ss’包含两个参数（$j和$j）；
#          threads->new(&ss,$j,$j);
#      }
#      foreach $thread(threads->list(threads::all))
#      {   if($thread->is_joinable())     #判断线程是否运行完成；
#          {   $thread->join();
#                #输出中间结果；
#              print scalar(threads->list()),"t$jt",localtime(time),"n";
#          }
#       }
# }




 
my @links=(
  ["http://www.foo.bar/rulez.data","rulez_data.txt"],
  ["http://new.host/more_data.doc","more_data.doc"],
);

 
# Max 30 processes for parallel download
my $pm = Parallel::ForkManager->new(30);
 
LINKS:
foreach my $linkarray (@links) {
  $pm->start and next LINKS; # do the fork
 
  my ($link, $fn) = @$linkarray;
  warn "Cannot get $fn from $link"
    if getstore($link, $fn) != RC_OK;
 
  $pm->finish; # do the exit in the child process
}
$pm->wait_all_children;