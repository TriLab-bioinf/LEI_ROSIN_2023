#!/Users/lorenziha/opt/anaconda3/envs/snp_caller/bin/perl
use strict;

my $usage = "$0 -b <bedtools gff3 file> -g <annotation gff3 file>\n\n";
my %arg = @ARGV;
die $usage unless $arg{-g} and $arg{-b};

# Load cluster data first

# PGA_scaffold_46__1_contigs__length_54420    GIGAdb  gene    13685   15383   .   +   .   ID="gene_16295"; Evidence_id="plodia_orig_RNA-NM_001160195.1_R2.p4";    93 
open(CLUST_GFF, "<$arg{-b}");
my (%GENES,%COORDS, %STRAND);
while(<CLUST_GFF>){
	chomp;
	my @x = split /\t/;
	my $gene_id = $1 if $x[8] =~ m/=.(gene_\d+)/;
	my ($lend, $rend, $feat, $cluster, $ctg, $strand) = ($x[3], $x[4], $x[2], $x[9], $x[0], $x[6]);  
	push @{$COORDS{$cluster}}, $rend, $lend;
	$STRAND{$cluster} = $strand;
	push @{$GENES{$cluster}}, $gene_id; 
	#print "Cluster number=$cluster\n";
}
close CLUST_GFF;

# sort new gene coords
my %genes;
my $counter;
foreach my $k (keys %COORDS){
	$counter++;
	@{$COORDS{$k}} = sort {$a <=> $b}  @{$COORDS{$k}};

        my $cluster_id = $k;
	my $lend = $COORDS{$k}[0];
	my $rend = $COORDS{$k}[-1];
        foreach my $g (@{$GENES{$k}}){
                #print "$g\t$coord\t$strand\t$cluster_id\n";
                $genes{$g}->{rend} = $rend;
                $genes{$g}->{lend} = $lend;
                $genes{$g}->{cluster_id} = $cluster_id; # print ">> $counter >>cluster_id=$genes{$g}->{cluster_id} -- g=$g -- cluster_id_var=$cluster_id -- k=$k\n";
		$genes{$g}->{strand} = $STRAND{$k};
        }
}


# process gff file
open(GFF, "<$arg{-g}");
my %cluster_count;
while(<GFF>){
	chomp;
	my @x = split /\t/;
	my $gene_id;
	if($x[8]=~m/=.(gene_\d+)/){
		$gene_id = $1;
	} else {
		# ID="mRNA_1"; Parent="gene_1" ;Evidence_id="plodia_orig_RNA-XM_038016986.1_R9.p1";
		print "$x[8] -> gene_id=$gene_id --\n";
		exit;
	}
	if($x[2] eq "gene"){
		# skip same cluster gene lines
		my $cluster = $genes{$gene_id}->{cluster_id};
		$cluster_count{$cluster}++;
		next if $cluster_count{$cluster} > 1;
		# replace gene_id with cluster id and rend and lend coords
		$x[3] = $genes{$gene_id}->{lend};
		$x[4] = $genes{$gene_id}->{rend};
		$x[8] =~ s/$gene_id/gene_N$genes{$gene_id}->{cluster_id}/;
	} elsif ($x[2] eq "mRNA"){
		# replace parent gene id
		$x[8] =~ s/$gene_id/gene_N$genes{$gene_id}->{cluster_id}/;
	} elsif ($x[2] eq "exon"){
		# replace gene_id
		$x[8] =~ s/$gene_id/gene_N$genes{$gene_id}->{cluster_id}/;
	} else {}

	print join("\t",@x)."\n";
}
close GFF;

