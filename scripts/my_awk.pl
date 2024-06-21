#!/usr/bin/perl

use Getopt::Long;
use strict; 
use warnings; 
# flag_v might have bug 
my $flag_v=0; 
GetOptions ('v' => \$flag_v);
my $list=$ARGV[0]; 
my $table=$ARGV[1];
my $col=$ARGV[2] // 0; # index  
my  $delim=$ARGV[3] // "\t";

my %hash; 
open (LIST,"<",$list) or die "Can't read $list!"; # Takes first column by tab 
while (<LIST>){
chomp; 
my @array=split("\t",$_);
my $item=$array[0]; 
$hash{$item}=1; 
} 
close LIST; 

open (TABLE,"<",$table) or die "Can't read $table!"; 
while (<TABLE>){
chomp; 
my @array=split($delim,$_); 
my $val=$array[$col];
if (exists $hash{$val}){
if ($flag_v !=1){
my $str=join($delim,@array); 
print "$str\n"; 
}
}
else{
if ($flag_v==1){
my $str=join($delim,@array); 
print "$str\n"; 
}
}
}
close TABLE; 
