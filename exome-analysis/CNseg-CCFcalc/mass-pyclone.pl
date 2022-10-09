#!/usr/bin/perl
# Analyze mutations with PyClone-VI
# Needs to be run on conda environment with PyClone-VI installed/loaded
use strict;
use warnings;
my $fh;
open ($fh,"<samples.txt") || die "c";
while (<$fh>){
	s/[\r\n]//g;
	my @lines = split ("\t",$_);
	my $f = $lines[0];
	my $pyf = $f.".pyclone.input.txt";
	system("pyclone-vi fit -i ./pyclone-input/$pyf -o ./pyclone-out/$f.h5 -c 40 -d beta-binomial -r 10");
	system("pyclone-vi write-results-file -i ./pyclone-out/$f.h5 -o ./pyclone-out/$f.pyclone.results.txt");

}
close $fh;
