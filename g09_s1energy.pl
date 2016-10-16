#!/usr/bin/perl
use strict;
use warnings;

my $searchphrase1 = "SCF Done";
my $searchphrase2 = "Optimization completed";
my $searchphrase3 = "Total Energy";
open(IN,"$ARGV[0]") or die "Cannot read $ARGV[0] file!";
my @lines;
while (<IN>) {
	if ((/$searchphrase1/) or (/$searchphrase2/) or (/$searchphrase3/)) {
	push @lines,$_;
	}
}	
my @states;
for my $i (0 .. $#lines) {
	if ($lines[$i] =~ /$searchphrase2/) {
		push @states,  "$lines[$_]\n"  for grep $_ >=0, $i-2 .. $i-1;	 
	}
}
close IN;
print "Ground state energy (eH): $states[0]";
print "Excited state energy (eH): $states[1]";
