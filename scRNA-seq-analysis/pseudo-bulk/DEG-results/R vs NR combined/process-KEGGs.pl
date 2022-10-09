#!/usr/bin/perl
use strict;
use warnings;
my $fh;

open ($fh,"<pre-select-KEGG.results.txt") || die "c";
<$fh>;
my %hash; #foreach combo of cell-group-pathway, saves the most significant one
while (<$fh>){
    s/[\r\n]//g;
    s/"//g;
    my $line = $_;
    my @lines = split ("\t",$line);
    my @blah = split (":",$lines[0]);
    my $cell = $blah[0];
    my @bb = split ("_",$blah[1]);
    my $group = $bb[0];
    my $dir = $bb[1];

    my $mult = 1;
    if ($dir eq "down") { $mult = -1; }

    my $nlp = $lines[10] * $mult; # signed qval
    #my $nlp = (log($lines[6])/log(10))*-1; # from p.adjust
    my $pathway = $lines[2];
    my $ratio = $lines[3];
    my @rr = split ("\/",$ratio);
    $ratio = $rr[0]/$rr[1] * $mult; # signed ratio
    my $index = "$cell:$group:$pathway";

    $hash{$index } = [$nlp,$ratio]; # will overwrite earlier entries for the sam cell-group-pathway;
    #  results file is already ordered by ascending significance
}
close $fh;

print "Cluster\tcellType\tgroup\tPathway\tnlp\tGeneRatio\n";
foreach my $k (sort keys %hash){
    my @blah = split (":",$k);
    print "$k\t$blah[0]\t$blah[1]\t$blah[2]\t$hash{$k}[0]\t$hash{$k}[1]\n";
}

