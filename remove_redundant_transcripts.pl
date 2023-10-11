#!/Users/lorenziha/opt/anaconda3/envs/snp_caller/bin/perl
use strict;

# output mRNA IDs to be removed

my %iso;
while(<>){
	chomp;
	my @x = split /\t/;
	next unless $x[2] eq 'exon';
	my $id = $1 if m/Parent=.(mRNA_\d+)/;
	$iso{$id} .= $x[0].'_'.$x[3].'_'.$x[4].'_';
}

my %nr_iso;
foreach my $k (keys %iso){
	$nr_iso{$iso{$k}}= $k;
}

foreach my $nk (keys %nr_iso){
	print  $nr_iso{$nk}."\n";
}
