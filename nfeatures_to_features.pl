$filename_nfeatures=@ARGV[0];
$filename_rscu=@ARGV[1];
#read rscu file
open(RSCU,$filename_rscu);
@rscu=<RSCU>;
$header_rscu=shift(@rscu);
$header_rscu=~s/\t/,/g;
$header_rscu=~s/\s+//g;
chop $header_rscu;

$rscu_values=shift(@rscu);
$rscu_values=~s/\t/,/g;
$rscu_values=~s/\s+//g;
chop $rscu_values;

close RSCU;
#read and process nfeature file
open(FH,$filename_nfeatures);
@file=<FH>;
close FH;
$header=shift(@file);

foreach $line(@file){
	$line=~s/,/kk,/g;	
	@line=split('kk',$line);
	$seq_id=shift(@line);
	$seq_id=~s/>//g;
	$new_header=$header;
	$new_header=~s/,/,$seq_id\_/g;
	$new_header=~s/Sequence_ID,//g;
	$all_header.=$new_header;
	#print @line;
	push(@array,@line);	
}
$all_header=~s/\n/,/g;
chop $all_header;
$array=join('',@array);
$array=~s/\n//g;
$array=~s/^,//g;
$feat_table="$all_header,$header_rscu\n"."$array,$rscu_values";
$outfile=@ARGV[2];
open(OUT,'>',$outfile);
print OUT $feat_table;
close OUT;
