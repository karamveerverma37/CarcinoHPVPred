#Created by Karamveer, DBT-BINC fellow, Pondicherry University on 05/01/2022 using Perl5 (v5.26.2)
#!/usr/bin/perl
#use Cwd qw(cwd);
#BEGIN {
#		$dir = cwd;
#	}
#use lib $dir;
#/home/karamveer/works/HPV/evolution/New/carc_vs_non_carc/DAMBE/modelling/CHCT';
#use RSCU_table;
if($#ARGV<1){
	die "Usage: nfeatures.pl input_file.fa output_file.tsv\n";
}
$filename=@ARGV[0];
$output_filename=@ARGV[1];
$path="/usr/local/bin/python3.9";
system("python3 Nfeature_DNA.py -i $filename -o $output_filename -ft ALL_COMP -path $path");
#system("python3 Nfeature_DNA.py -i HPV18_out.fa -o nfeature_out.tsv -ft ALL_COMP -path /usr/local/bin/python3.9");
#"E2_PDNC_TC" "E2_CDK_TA" "E2_RDK_TA" "E6_ENT_NL_A" "E1_RDK_AG" "E2_PDNC_AA" "E1_PKNC_AAA" "E1_PDNC_AA" "E2_PKNC_AAA" "E2_PKNC_ATC" "E2_RDK_AG" "E1_DDON_A" "E1_A3" "L2_PDNC_AA" "E2_PKNC_TCA" "E2_Entropy" "E1_PKNC_CTT" "E2_CDK_CT" "E2_ENT_NL_A" "E6_PKNC_GAA" "L2_PKNC_AAA" "L2_PDNC_GC" "E1_PDNC_GA" "L1_PDNC_GC" "E1_PKNC_GGA" "E2_PDNC_GA" "E6_PKNC_AAA" "E1_CDK_CT" "E6_DDON_A" "L2_CDK_GA" "L1_PKNC_GAT" "E2_PDNC_AT" "E2_PDNC_TG" "E2_PKNC_GAA" "E2_PKNC_GTC" "E2_CDK_AG" "E2_PKNC_GGA" "E6_PDNC_AA" "E2_RDK_CA" "E1_PKNC_ACC")

=pod
%sequences;
open( FH, '<', $filename ) or die $!;
while($line=<FH>){
      chomp $line;
      if($line=~/^(>.*)$/){
           $id = $1;
      }
      elsif ($line !~ /^\s*$/){
           $sequences{$id} .= $line;
      }
}
foreach $key (keys(%sequences)){
	$seq=$sequences{$key};
	$fas_seq="$key\n$seq";
	print "$fas_seq\n";
	#system("python3 Nfeature_DNA.py -i  -o COMP -ft ALL_COMP")
}
=cut
