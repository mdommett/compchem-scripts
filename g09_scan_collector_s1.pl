#!/usr/bin/perl
use strict;
use warnings;

my $searchphrase1 = "SCF Done";
my $searchphrase2 = "Optimization completed";
my $searchphrase3 = "Total Energy";
open(IN,"$ARGV[0]") or die "Cannot read $ARGV[0] file!";
open(OUT, ">tmp.txt") or die "Cant read destination file: $!\n";
my @lines;
while (<IN>) {
	if ((/$searchphrase1/) or (/$searchphrase2/) or (/$searchphrase3/)) {
	push @lines,$_;
	}
}	


for my $i (0 .. $#lines) {
	if ($lines[$i] =~ /$searchphrase2/) {
		printf OUT  "$lines[$_]\n"  for grep $_ >=0, $i-2 .. $i-1;	 
	}
}
close IN;

open(FINAL, ">$ARGV[0].energies.txt") or die "Cant open final: $!\n";

my @g;
my @S0;
my @S1;
my $eVconv = 27.2114;
open(OUT, "<tmp.txt") or die "Cant read destination file: $!\n";

printf FINAL "S0 S1 \n";

while (<OUT>) {
	if (/$searchphrase1/) {
		chomp;
		$_=~s/^\s*//;
		$_=~s/^\s*$//;
		@g=split(/\s+/,$_);
		push @S0, $g[4];
#		printf FINAL "S0 energy eV: % .8f \n", $S0*$eVconv;
	}
	
	if (/$searchphrase3/) {
		chomp;
                $_=~s/^\s*//;
                $_=~s/^\s*$//;
                @g=split(/\s+/,$_);
                push @S1, $g[4];
 #               printf FINAL "S1 energy is: % .8f \n\n", $S1*$eVconv;
	}
}

for (my $i=0; $i <= $#S0; $i++) {
   printf FINAL "% .8f % .8f \n", $S0[$i]*$eVconv-$S0[0]*$eVconv, $S1[$i]*$eVconv-$S0[0]*$eVconv
}

 
				

			
			

		
