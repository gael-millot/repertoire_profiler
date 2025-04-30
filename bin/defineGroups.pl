#!/usr/bin/perl -w

use strict;
use warnings;

sub groups{
    my ($igcall) = @_;

    my %groups;
    
    my @cols = split(/,/,$igcall);
    for my $c (@cols){
	$c=~s/\*.*$//;
    $c =~ s/\//_/g; # replace any "/" characters by "_" because files will be named with this
	$groups{$c}=1;
    }
    return %groups;
}


<>;
while(<>){
    chomp;
    my @cols = split(/\t/);

    my $name= $cols[0];
    my $ig1 = $cols[4];
    my $ig2 = $cols[6];

    my %groups1=groups($ig1);
    my %groups2=groups($ig2);

    for my $g1 (keys(%groups1)){
	for my $g2 (keys(%groups2)){
	    print "$name\t--$g1--$g2--\n";
	}
    }
}
