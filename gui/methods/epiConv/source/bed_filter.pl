#!/usr/bin/perl -w

##filter bed file given a list of cell barcodes
##$ARGV[0]: cell barcode file
##$ARGV[1]: bed file

use strict;

open IN, "<$ARGV[0]";

my %barcode;
my $i=1;
while(<IN>){
	chomp;
	my @ident=split /\t/, $_;
	$barcode{$ident[0]}=$i;
	$i++;
}
close IN;

open IN, "<$ARGV[1]";
while(<IN>){
	chomp;
	my @record=split /\t/, $_;
	if(exists $barcode{$record[3]}){
		print "$record[0]\t$record[1]\t$record[2]\t$barcode{$record[3]}\n";
	}
}
close IN;
