#!/usr/bin/perl
use strict;
use warnings;
my $fh;

my %hash;
open ($fh,"<valero-ng-data.annot.txt") || die "c";
<$fh>;
while (<$fh>){
    s/[\r\n]//g;
    my $line = $_;
    my @lines = split ("\t",$line);
    $hash{$lines[0]} = [$lines[19],$lines[21]];
}

open ($fh,"<valero-nc-data.txt") || die "c";
my $head = <$fh>;
$head =~ s/[\r\n]//g;
print "$head\tMMR_canon\tARID1A_mut\n";
while (<$fh>){
    s/[\r\n]//g;
    my $line = $_;
    my @lines = split ("\t",$line);
    if (defined $hash{$lines[0]}){
        print "$line\t$hash{$lines[0]}[0]\t$hash{$lines[0]}[1]\n";
    }
}
close $fh;