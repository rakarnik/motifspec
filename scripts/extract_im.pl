#!/usr/bin/perl

use strict;

my $input = $ARGV[0];
my $motif = $ARGV[1];

open INPUT, "< $input"
    or die "Unable to open ACE file '$input'\n";
my ($line, @fields, $seqs, $seqnum, $uniq, $in_site_list);
$seqnum = 0;
while (<INPUT>) {
    chomp;
		next if /^$/;

		if(/^Motif $motif/) {
			#print "In motif $motif\n";
			$in_site_list = 1;
			next;
    }
    if(/^\*/ && $in_site_list) {
			#print "Leaving motif $motif\n";
			$in_site_list = 0;
			last;
    }
    
    if(/^\#/) {
			@fields =split /\s+/;
			$seqs->{"$seqnum"} = $fields[1];
			$seqnum++;
    }
    if($in_site_list) {
			@fields = split /\s+/;
			$uniq->{$seqs->{@fields[1]}} = 1;
    }
} 
close INPUT;

my $gene;
foreach $gene(keys %$uniq) {
    print "$gene\n";
}
