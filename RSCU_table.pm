#!/usr/bin/perl
package RSCU_table;
my $LEVEL=1;
sub table{
	$filename=$_[0];
	open(FH,$filename);
	@file=<FH>;
	$file1=join('',@file);
	@file2=split('>',$file1);
	shift @file2;
	foreach $line(@file2){ 
		$gene_id = substr($line,0,2,'');
		@rscu1= RSCU($line);
		push(@array1,"$gene_id\n@rscu1\ncutpart");
#		@rscu1=();
	} 
	return @array1;
}

#@rscu2=();
#@output=();
sub RSCU{
	my $seq = $_[0];
	$seq=~s/\s+//g;
	my @output=();
	#print "$seq\n";
	$seq =~ tr/actguU/ACTGTT/;
	#print "$seq\n";
	my %codigenetic = ("TTT","F","TCT","S","TAT","Y","TGT","C","TTC","F","TCC","S","TAC","Y","TGC","C","TTA","L","TCA","S","TAA",".","TGA",".","TTG","L","TCG","S","TAG",".","TGG","W","CTT","L","CCT","P","CAT","H","CGT","R","CTC","L","CCC","P","CAC","H","CGC","R","CTA","L","CCA","P","CAA","Q","CGA","R","CTG","L","CCG","P","CAG","Q","CGG","R","ATT","I","ACT","T","AAT","N","AGT","S","ATC","I","ACC","T","AAC","N","AGC","S","ATA","I","ACA","T","AAA","K","AGA","R","ATG","M","ACG","T","AAG","K","AGG","R","GTT","V","GCT","A","GAT","D","GGT","G","GTC","V","GCC","A","GAC","D","GGC","G","GTA","V","GCA","A","GAA","E","GGA","G","GTG","V","GCG","A","GAG","E","GGG","G");

	#print "$codigenetic{'TTT'}";
	my %Codo = ('TTT',0, 'TTC',0, 'TTA',0, 'TTG',0, 'CTT',0, 'CTC',0, 'CTA',0, 'CTG',0, 'ATT',0, 'ATC',0, 'ATA',0, 'ATG',0, 'GTT',0, 'GTC',0, 'GTA',0, 'GTG',0, 'TCT',0, 'TCC',0, 'TCA',0,'TCG',0, 'CCT',0, 'CCC',0, 'CCA',0, 'CCG',0, 'ACT',0, 'ACC',0, 'ACA',0, 'ACG',0, 'GCT',0, 'GCC',0,'GCA',0, 'GCG',0, 'TAT',0, 'TAC',0, 'TAA',0, 'TAG',0, 'CAT',0, 'CAC',0, 'CAA',0, 'CAG',0, 'AAT',0,'AAC',0, 'AAA',0, 'AAG',0, 'GAT',0, 'GAC',0, 'GAA',0, 'GAG',0, 'TGT',0, 'TGC',0, 'TGA',0, 'TGG',0,'CGT',0, 'CGC',0, 'CGA',0, 'CGG',0, 'AGT',0, 'AGC',0, 'AGA',0, 'AGG',0, 'GGT',0, 'GGC',0, 'GGA',0,'GGG',0);

	#Codon usage calculation
	my $seq2 = $seq;
	$seq2 =~ s/(...)/$1 /g;
	#print "$seq2\n";
	my %totalCodo =();
	my %totalAA =();
	my %CodonsPerAA;
	my %rscu;
	my $AA='';
	my $totalCodons = 0;
	foreach my $codo(sort keys(%Codo)) {
		my $num = $seq2 =~ s/$codo /$codo /g;
		#print "$codo\t$num\n";
		$num += 0;
		$totalCodo{$codo} = $num;
		$AA = $codigenetic{$codo};
		$totalAA{$AA} += $num;
		$totalCodons += $num;
	}

	#print %totalCodo;
	#print %totalAA,"\n";
	#print "$totalCodons\n";
	#RSCU calculation
	#%rscu = ();



	foreach my $codon(keys(%codigenetic)) {
		if ($codon ne ""){
			$CodonsPerAA{$codigenetic{$codon}}++;
		}
	}
	foreach $AA(keys(%CodonsPerAA)) {
		push @CodonsPerAA_aux, $AA;
		push @CodonsPerAA_aux, $CodonsPerAA{$AA};
	}


	foreach $codo(sort keys(%Codo)) {
		$AA = $codigenetic{$codo};
		if ($totalAA{$AA} == 0) {
			$rscu{$codo} = 0;
		}else {
			$rscu{$codo} = ($totalCodo{$codo}  / $totalAA{$AA} ) * $CodonsPerAA{$AA};
		}
	}
	foreach my $key(keys %rscu){
		if(($key ne 'TGA') && ($key ne 'TAA') && ($key ne 'TAG')){
			push(@output,"$key\t$rscu{$key}\n");
		}
	}
	return @output;
}
1;
