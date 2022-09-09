#Created by Karamveer, DBT-BINC fellow, Pondicherry University on 05/01/2022 using Perl5 (v5.26.2)
#!/usr/bin/perl
use Cwd qw(cwd);
BEGIN {
		$dir = cwd;
	}
use lib $dir;
#/home/karamveer/works/HPV/evolution/New/carc_vs_non_carc/DAMBE/modelling/CHCT';
use RSCU_table;
$filename=@ARGV[0];
if($#ARGV<1){
	die "Usage: perl Codon_analysis.pl input_file.fa output_file.tsv\n";
}
@table=RSCU_table::table($filename);
$rscu_joined = join('',@table);
@each_gene_rscu=split('\ncutpart',$rscu_joined);
#print "***$each_gene_rscu[1]***\n";

#open(FH,$filename);
#@file=<FH>;

#@output=("Gene\t","A3\t","T3\t","G3\t","C3\n");
foreach $gene(@each_gene_rscu){
	#chop $gene;
	@keys=split('\n', $gene);
	#print "$keys[1]\n";	
	$gene_ids=shift(@keys);
	#print "$gene_ids\n";
	$value_A=0; $value_T=0; $value_G=0; $value_C=0;
	%hash_table=();
	foreach $key(@keys){
		#print "***$key***\n";
		@codon_usage=split('\t',$key);
		$codons=$codon_usage[0];
		$rscu_value=$codon_usage[1];
		#print "$codons\t$rscu_value\n";
		if($rscu_value>=1.6){
			$hash_table{$codons}=$rscu_value;
		}
	}
	foreach my $new_codon (keys %hash_table){
		my $value=$hash_table{$new_codon};
		#print "$new_codon\t$value\n";
		$new_codon=~s/\s+//g;
		@codon=split('',$new_codon);
		#print"***$codon[2]***\n";
		if($codon[2] eq 'A'){
			$value_A += $value;
		}
		elsif($codon[2] eq 'T'){
                	$value_T += $value;
                }
		elsif($codon[2] eq 'G'){
                        $value_G+=$value;
                }
		elsif($codon[2] eq 'C'){
                        $value_C+=$value;
                }

	}
	push(@header, "$gene_ids\_A3\t$gene_ids\_T3\t$gene_ids\_G3\t$gene_ids\_C3\t");
	push(@output,"$value_A\t$value_T\t$value_G\t$value_C\t");
}
push(@header,"\n@output");
#print "@header\n";
chomp(@header);
open(OUT, '>', $ARGV[1]) or die "provide a output filename.\n";
print OUT "@header\n";
#print "Output RSCU table file $ARGV[1] created successfully.\n\n";
