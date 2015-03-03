#!/usr/bin/perl
use strict;

my ($line, @f);

printf "%-20s %5s %6s %6s %6s %6s %8s %-35s %7s %5s %8s\n", 
				"ace", "mot", "tot", "ssp", "gtseq", "hits", "score", "cons", "seqc", "sspc", "iter";
my ($file, $filenopath, $tot, $mot, $seen); 
my ($ssp, $gtseq, $hits, $score, $seqcut, $sspcut, $iter, $consensus);
foreach $file (@ARGV) {
	open ACE, "<$file";
	$filenopath = $file;
	$filenopath =~ s/.*\///;
	$tot = 0;
	while (chomp($line=<ACE>)) {
		my @cons;
		if ($line =~ /^\#/) {
			$tot++;
		}
	  if ($line =~ /Motif/) {         # Start reading motif
			@f = split /\s+/, $line;
			$mot = $f[1];
			$seen = {};
			$hits = 0;
			chomp($line = <ACE>);
			while (!($line =~ /\*/)) {    # Read until ******...
		    @f = split /\s+/, $line;
		    push @cons, $f[0];
				if(! defined $seen->{"$f[1]"}) {
					$hits++;
					$seen->{"$f[1]"} = 1;
				}
		    chomp($line = <ACE>);
			}                              # Done reading motif here
				
			chomp($line = <ACE>);          # Read score
			@f=split /:/, $line;
			$score = $f[1];

			chomp($line = <ACE>);          # Read number of sequences above sequence threshold
			@f=split /:/, $line;
			$gtseq = $f[1];
			
			chomp($line = <ACE>);          # Read size of search space
			@f=split /:/, $line;
			$ssp = $f[1];

			chomp($line = <ACE>);          # Read sequence cutoff
			@f=split /:/, $line;
			$seqcut = $f[1];

			chomp($line = <ACE>);          # Read search space cutoff
			@f=split /:/, $line;
			$sspcut = $f[1];

			chomp($line = <ACE>);          # Read iteration found
			@f=split /\s+/, $line;
			$iter = $f[2];
					
			chomp($line = <ACE>);		       # Eat empty line
	
			# Calculate consensus
			$consensus = "";
			my (%ct, $n, $l);
			foreach my $i (0..length($cons[0])-1) {
				$ct{"A"}=0;
				$ct{"C"}=0;
				$ct{"G"}=0;
				$ct{"T"}=0;
				foreach my $s (@cons) {
					$ct{substr($s, $i, 1)}++;
				}
				#print "a=", $ct{"A"},"\tc=", $ct{"C"}, "\tg=", $ct{G}, "\tt=", $ct{"T"}, "\n";

				my $max = max($ct{"A"}, $ct{"C"}, $ct{"G"}, $ct{"T"});
				$n=0;
				$ct{"A"}=$ct{"A"}/$max;
				$ct{"C"}=$ct{"C"}/$max;
				$ct{"G"}=$ct{"G"}/$max;
				$ct{"T"}=$ct{"T"}/$max;
				if ($ct{"A"}>0.5) {$n++;}
				if ($ct{"C"}>0.5) {$n++;}
				if ($ct{"G"}>0.5) {$n++;}
				if ($ct{"T"}>0.5) {$n++;}
				if ($n==1) {
					if ($ct{"A"}>0.5) {$l="A";}
					if ($ct{"C"}>0.5) {$l="C";}
					if ($ct{"G"}>0.5) {$l="G";}
					if ($ct{"T"}>0.5) {$l="T";}
				} elsif ($n==2) {
					if (($ct{"A"}>0.5)&&($ct{"C"}>0.5)) {$l="M";}
					if (($ct{"G"}>0.5)&&($ct{"T"}>0.5)) {$l="K";}
					if (($ct{"A"}>0.5)&&($ct{"G"}>0.5)) {$l="R";}
					if (($ct{"T"}>0.5)&&($ct{"C"}>0.5)) {$l="Y";}
					if (($ct{"A"}>0.5)&&($ct{"T"}>0.5)) {$l="W";}
					if (($ct{"G"}>0.5)&&($ct{"C"}>0.5)) {$l="S";}
				}	else {
					$l="-";
				}
				$consensus = $consensus . $l;
			}

			printf "%-20s %5d %6d %6d %6d %6d %8.2f %-35s %7.4f %5.2f %8s\n",
						 $filenopath, $mot, $tot, $ssp, $gtseq, $hits, $score, 
						 $consensus, $seqcut, $sspcut, $iter;
		}
	}
}
	
# Subroutine max()
sub max { 
  my ($a,$b,$c,$d)=@_;
  my $x=$a;
  if ($b>$x) {$x=$b;}
  if ($c>$x) {$x=$c;}
  if ($d>$x) {$x=$d;}
  $x;
}
