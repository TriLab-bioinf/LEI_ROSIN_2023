### **Generation of annotation gff file for new *P. interpunctella* assembly**

Annotation source: <http://gigadb.org/dataset/view/id/102231/File_sort/type_id>

*Protocol:*

*P. interpunctella* predicted transcript sequences from GIGAdb were generated from *plo_final_annotation.gff* using the following commands:

```         
sed -e 's/prediction/mRNA/' plo_final_annotation.gff > plo_final_annotation_mRNA.gff

gffread plo_final_annotation_mRNA.gff -g PinterpunctellaAssembly.asm.p_ctg_editedv6.fa -w plo_final_annotation.transcripts.fasta -A
```

Note: *gffread* utility was downloaded from <https://github.com/gpertea/gffread>

Then transcripts sequences were mapped to the new *P. interpunctella* assembly with *minimap2* (v2.26):

```         
minimap2 -ax splice --cs Plodia_genome_Scully_2022-edit.fa plo_final_annotation.transcripts.fasta > minimap2.alignments.plo_final.tmp.sam
```

The resulting sam file was sorted by read coordinate:

```         
samtools sort -o minimap2.alignments.plo_final.bam minimap2.alignments.plo_final.tmp.sam
```

The sorted bam file was processed with *sam_to_gtf.pl* to generate a temporary gtf file:

```         
./sam_to_gtf.pl minimap2.alignments.plo_final.bam > minimap2.alignments.plo_final.gtf
```

Finally, *minimap2.alignments.gtf* was processed by *get_genes_from_gtf.pl* script to generate the final gff file with the annotation:

```         
cat minimap2.alignments.plo_final.gtf|./get_genes_from_gtf.pl > Pinterpunctella_LEAH.plo_final.gff
```

**New gene merge protocol by percent overlap (\>=50% of smallest gene)**

Identification of GFF genes that overlap at least 50% of their length:

```         
./merge_genes_from_gff.py -g Pinterpunctella_LEAH.plo_final.gff -o 0.5 > tmp.500bp_NEW.bed
```

Then, overlapping genes were flagged as isoforms as described below:

```         
./merge_overlapping_genes.pl -b tmp.500bp_NEW.bed -g Pinterpunctella_LEAH.plo_final.gff >  Pinterpunctella_LEAH.plo_final.clustered_NEW.gff.tmp
```

The resulting temp gff file was then sorted by gene and mRNA position:

```         
egrep '\tgene\t' Pinterpunctella_LEAH.plo_final.clustered_NEW.gff.tmp |cut -f 2 -d '"' > cluster_NEW.ids

for i in `cat cluster_NEW.ids`
     do echo $i
     grep -w $i Pinterpunctella_LEAH.plo_final.clustered_NEW.gff.tmp >> sorted_NEW.gff
     done
```

The sorted gff file was afterwards processed to eliminate redundant isoforms that shared exactly the same exonic structure:

```         
cat sorted_NEW.gff| ./remove_redundant_transcripts.pl > mRNA_NEW.NR.ids

egrep -wf mRNA_NEW.NR.ids sorted_NEW.gff > sorted_NEW.NR.gff
```
### **Identification of syntenic blocks between new *P. interpunctella* (GCF_027563975.1) and *Bombyx mori* (GCF_014905235.1) assemblies** 

*P. interpunctella* and *Bombyx mori* assemblies were aligned in protein space with the *mummer* tool *promer* and the resulting delta file was formatted as a hit coordinate file.  
```
promer -p Pinter_vs_Bmori_promer GCF_014905235.1_Bmori_2016v1.0_genomic.fna GCF_027563975.1_ilPloInte3.1_genomic.fna

show-coords -clrT Pinter_vs_Bmori_promer.delta > Pinter_vs_Bmori_promer.coords
```

Next, the *Pinter_vs_Bmori_promer.coords* file was parsed as a matchList file, to be used as input for the tool *DAGchainer* (https://github.com/kullrich/dagchainer).   

```
cat Pinter_vs_Bmori_promer.coords|grep '^[0-9]' >  Pinter_vs_Bmori_promer.NO_HEADER.coords
 
cat Pinter_vs_Bmori_promer.NO_HEADER.coords | \
     perl -lane '$gene_1="$F[15]_$F[0]_$F[1]"; $gene_2="$F[16]_$F[2]_$F[3]" ;($int,undef)=split(/\./,$F[6]);print "$F[15]\t$gene_1\t$F[0]\t$F[1]\t$F[16]\t$gene_2\t$F[2]\t$F[3]\t1e-$int\t$F[7]"' > Pinter_vs_Bmori_promer.matchList
```

File *Pinter_vs_Bmori_promer.matchList* was then processed with the *DAGchainer* sofwtare to identify syntenic blocks between the two genomes
```
./accessory_scripts/filter_repetitive_matches.pl 5 < Pinter_vs_Bmori_promer.matchList > Pinter_vs_Bmori_promer.matchList.filtered

run_DAG_chainer.pl -i Pinter_vs_Bmori_promer.matchList.filtered -D 50000 -A 6 > Pinter_vs_Bmori_promer.matchList.filtered.synteny
```
