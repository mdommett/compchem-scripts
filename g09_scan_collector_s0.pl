#!/usr/bin/perl
use strict;
use warnings;

my $searchphrase1 = "SCF Done";
my $searchphrase2 = "Optimization complete";

open(IN,"$ARGV[0]") or die "Cannot read $ARGV[0] file!";
open(OUT, ">tmp.txt") or die "Cant read destination file: $!\n";
my @lines;
while (<IN>) {
	if ((/$searchphrase1/) or (/$searchphrase2/)) {
	push @lines,$_;
	}
}	


for my $i (0 .. $#lines) {
	if ($lines[$i] =~ /$searchphrase2/) {
		printf OUT  "$lines[$_]\n"  for grep $_ >=0, $i-1;
	}
}
close IN;

open(FINAL, ">$ARGV[0].energies.txt") or die "Cant open final: $!\n";

my @g;
my @S0;
my @S1;
my $eVconv = 27.2114;
open(OUT, "<tmp.txt") or die "Cant read destination file: $!\n";


while (<OUT>) {
	if (/$searchphrase1/) {
		chomp;
		$_=~s/^\s*//;
		$_=~s/^\s*$//;
		@g=split(/\s+/,$_);
		push @S0, $g[4];
	}
	
}

for (my $i=0; $i <= $#S0; $i++) {
   printf FINAL "% .9f \n", $S0[$i]
}

 
				

			
			

		
