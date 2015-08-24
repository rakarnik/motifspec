#!/usr/bin/perl

use strict;

my $ace = shift;
my $gcback = shift;
$gcback = 0.5 unless $gcback;

my $atback = (1 - $gcback)/2;
$gcback /= 2;

my $nt = {
	"A" => 0,
	"C" => 1,
	"G" => 2,
	"T" => 3
};

print "MEME version 4\n\n";
print "ALPHABET= ACGT\n\n";
print "strands: + -\n\n";
print "Background letter frequencies\n";
printf "A %0.2f C %0.2f G %0.2f T %0.2f\n\n", $atback, $gcback, $gcback, $atback;

open ACEIN, "< $ace";
my ($row, @f, @seq, $nseq, @mat, $col, $nmot);
$nmot = 1;
while(<ACEIN>) {
	next unless /Motif/;
	@seq = ();
	$nseq = 0;
	while(<ACEIN>) {
		chomp;
		last if /\*/;
		@f = split /\s+/;
		push @seq, $f[0];
		$nseq++;
	}
	$row = <ACEIN>;
	chomp $row;
	@f = split /\s+/, $row;
	print "MOTIF M$nmot $f[1]\n";
	@mat = ();
	$col = 0;
	foreach(@seq) {
		@f = split //;
		$col = 0;
		foreach(@f) {
			$mat[$col * 4 + $nt->{"$_"}]++;
			$col++;
		}
	}	
	print "letter-probability matrix: alength= 4 w= $col nsites= $nseq E= 0.0001\n";
	for(my $i = 0; $i < $col; $i++) {
		for(my $j = 0; $j < 4; $j++) {
			print " " unless $j == 0;
			printf "%1.6f",  $mat[$i * 4 + $j]/$nseq;
		}
		print "\n";
	}
	print "\n\n";	
	$nmot++;
}
close ACEIN;
