#!/home/lorenziha/data/miniconda3/bin/perl
use strict;

my %genes;
my %strand;
my %chr;
my %cds;
my  @my_keys;
my %nr;

while(<>){
    my @x = split /\t/;
    my $id = $1 if  $x[8] =~ m/^ID=(.+):?\;/;
    # add cds coordinates to %genes
    push @{$genes{$id}}, $x[3], $x[4];
    $strand{$id} = $x[6];
    $chr{$id} = $x[0];

    # cds
    my @coord = ($x[3], $x[4]);
    push @{$cds{$id}}, \@coord;  

    unless (exists($nr{$id})){
        push @my_keys, $id;
        $nr{$id}++;
    }
}

my $gene_id;
my $counter;
foreach my  $id (@my_keys){
    # print "MY ID=$id\n";
    $counter++;    
    $gene_id = "gene_".$counter;     
    my @sorted = sort {$a <=> $b} @{$genes{$id}};
    my $gene_start =  $sorted[0];
    my $gene_end = $sorted[-1];

    # Print gene line
    print join("\t", ($chr{$id},"GIGAdb","gene", $gene_start, $gene_end, ".", $strand{$id}), ".", "ID=\"$gene_id\"; Evidence_id=\"$id\";\n");

    # Print transcript line (mRNA)
    my $transcript_id = "mRNA_".$counter; 
    print join("\t", ($chr{$id},"GIGAdb","mRNA", $gene_start, $gene_end, ".", $strand{$id}), ".", "ID=\"$transcript_id\"; Parent=\"$gene_id\" ;Evidence_id=\"$id\";\n");

    # Print CDS line
    my @sorted_cds;
    my $cdsid;
    if($strand{$id} eq "+"){
        $cdsid = 0;     
        @sorted_cds =  @{$cds{$id}}; # keep the order
    } else {
        $cdsid = scalar(@{$cds{$id}}) + 1;    
        @sorted_cds =  reverse(@{$cds{$id}}); # reverse the order
    } 
    

    foreach my $cds_coords (@{$cds{$id}}){
        $cdsid = $strand{$id} eq "+"? $cdsid + 1 : $cdsid - 1 ;    
        my ($lend, $rend) = @{$cds_coords};
        my $exonid = "exon_$counter.$cdsid"; 
        print join("\t", ($chr{$id},"GIGAdb","exon", $lend, $rend, ".", $strand{$id}), ".", "ID=\"$exonid\"; Exon_number=\"$cdsid\"; gene_id=\"$gene_id\";Parent=\"$transcript_id\"\n");
    } 

}


 
