#!/usr/bin/perl
use Cwd qw(cwd);
BEGIN {
		$w_dir = cwd;
	}
use lib $w_dir;
#/home/karamveer/works/HPV/evolution/New/carc_vs_non_carc/DAMBE/modelling/CHCT';
use RSCU_table;

$genome_file_name=@ARGV[0];
$output_file_name=@ARGV[1];
if($#ARGV<1){
	die "Usage: perl Genome_to_gene.pl input_file.fasta output.fa\n";
}
#print $genome_file_name;
@dir=split('\.',$genome_file_name);
$dir= $dir[0];
if (-d "$dir") {
	system("rm -r $dir")
}
system("prokka/bin/prokka --kingdom Viruses  $genome_file_name --outdir $dir --prefix $dir --quiet");

$predicted_gene_file="$dir/$dir.ffn";
#print $predicted_gene_file;
open(FILE,$predicted_gene_file);
@file=<FILE>;
$file1=join('',@file);
$file1=~s/>/cutpart>/g;
$file1=~s/>.*[Pp]rotein E1/>E1/g;
if($file1=~">.*[Pp]rotein E2"){
	$file1=~s/>.*[Pp]rotein E2/>E2/g;
	}
	else
	{
	$file1=~s/>.*hypothetical protein/>E2/g;
	}

#$file1=~s/>.*[Pp]rotein E2/>E2/g;
if($file1=~">.*[Pp]rotein E6"){
	$file1=~s/>.*[Pp]rotein E6/>E6/g;
	}
	else
	{
	$file1=~s/>.*hypothetical protein/>E6/g;
	}
#$file1=~s/>.*[Pp]rotein E6/>E6/g;
$file1=~s/>.*[Pp]rotein E7/>E7/g;
$file1=~s/>.*[Pp]rotein L1/>L1/g;
if($file1=~">.*[Pp]rotein L2"){
	$file1=~s/>.*[Pp]rotein L2/>L2/g;
	}
	else
	{
	$file1=~s/>.*hypothetical protein/>L2/g;
	}

@file1=split('cutpart',$file1);
#print @file1;
foreach $line(@file1){
	if(($line=~">E1\n")||($line=~">E2\n")||($line=~">E6\n")||($line=~">E7\n")||($line=~">L1\n")||($line=~">L2\n")){	
		push(@new_file,$line);
	}
}

open(OUT, '>',$output_file_name) or die $!;
print OUT @new_file;
close OUT;
