#! /usr/bin/env perl 
use strict;
use warnings;
use Data::Dumper;

=head
-----------------------------------------------------

This script formats TopHat2 outputs into files 
than can be directly taken in by 'PhyExp'

Input: a list of TopHat2 output UniqueReads files
and a file contains orthlog gene ids

Output: a merged file of all the reads count data
and a file contains orthlog gene ids with corresponding gene lengths 

Author: Hang Ruan <hang.ruan@hotmail.com>

-----------------------------------------------------
=cut


my $lst = shift || die "$0 <tophat2out.lst> <orthlog.txt> 1>readsCount.txt 2>geneInfo.txt\n";
my $orthlog = shift || die "$0 <tophat2out.lst> <orthlog.txt> 1>readsCount.txt 2>geneInfo.txt\n";

my %readsCount;

open LST, $lst || die "Can not open $lst\n";

my @bio_rep;

my %gene_length;

while (<LST>) {
	
	chomp;
	
	open FP, $_ || die "Can not open $_\n";
	
	my $taxa_id = (split /\_/, $_)[0];

	my $head = <FP>;
	
	$head =~ s/\"//g;

	my @tab = split /\s+/, $head;
	my @bio_rep_id;	

	for (my $i=6; $i<= $#tab; $i++) {
		push @bio_rep_id, $tab[$i];
	}	

	push @bio_rep, @bio_rep_id;
	
	while (<FP>) {
		s/\"//g;
		my @tab = split;
	
		my $gene_id = $tab[0];
		my $exon_length = $tab[5];

		$gene_length{$gene_id} = $exon_length;
	
		for (my $i=6;$i<= $#tab;$i++) {
			$tab[$i] *= $exon_length;
			$readsCount{$taxa_id}{$gene_id}{$bio_rep_id[$i-6]} = int($tab[$i]+0.5);
		}	
	}

	close FP;	
	
}
	
close LST;


#print Dumper(\%readsCount);
1;

print "homoSapienGeneId";

foreach (@bio_rep) {
	print "\t$_";
}

print "\n";

1;

open ORTH, $orthlog || die "Can not open $orthlog\n";

my $head = <ORTH>;

my @taxa_id = split /\s+/, $head;

#print STDERR "Human_Gene_Id";

for (my $i=0;$i<$#taxa_id;$i++) {
	print STDERR "$taxa_id[$i]\t";
}

print STDERR "$taxa_id[-1]\n";

#my %bio_rep = map {$_ =>1} @bio_rep;

while (<ORTH>) {

	my @orth_gene_id = split;

	print "$orth_gene_id[0]";	

	#print STDERR "$orth_gene_id[0]";
	
	for (my $i=0;$i<=$#orth_gene_id;$i++) {

		if ($i == $#orth_gene_id) {
			print STDERR "$orth_gene_id[$i]:$gene_length{$orth_gene_id[$i]}";
		} else {
			print STDERR "$orth_gene_id[$i]:$gene_length{$orth_gene_id[$i]}\t";
		}	
	
		foreach my $rep(@bio_rep) {
			next unless ($rep =~ /$taxa_id[$i]/);
			if (exists $readsCount{$taxa_id[$i]}{$orth_gene_id[$i]}{$rep}) { 
				print "\t$readsCount{$taxa_id[$i]}{$orth_gene_id[$i]}{$rep}";
			} else {
			}
		}
			
	}
	
	print STDERR "\n";
	
	print "\n";		
}

close ORTH;


