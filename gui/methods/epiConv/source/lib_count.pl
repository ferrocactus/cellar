#!/usr/bin/perl -w
use strict;

open IN, "<$ARGV[1]";
my %res;
my @barcode;
while(<IN>){
	chomp;
	$res{$_}=0;
	push @barcode,$_;
}
close IN;

open IN, "<$ARGV[0]";
while(<IN>){
	chomp;
	my @record=split /\t/,$_;
	my $barcode=$record[3];
	$res{$barcode}++;
}
close IN;

while(my $key=shift @barcode){
	print "$key\t$res{$key}\n";
}

