if($#ARGV<2){
	die "Usage: Tool_all.pl HPV.fasta model input_feature_type";
}
$input_genome=@ARGV[0];
$model=@ARGV[1];
$input_type=@ARGV[2];
if($input_type=~m/full_genome/){
	@input_genome1=split('\.',$input_genome);
	$prefix=$input_genome1[0];
	$final_features= "$prefix".'_features.csv';
	$out_genes="$prefix".'_genes.fa';
	$tmp_out_features='HPV_features.tsv';
	$tmp_out_rscu='HPV_rscu.tsv';
	system("perl genome_to_gene.pl $input_genome $out_genes");
	system("rm -f $input_genome");
	system("perl nfeatures.pl $out_genes $tmp_out_features");
	system("perl codon_analysis.pl $out_genes $tmp_out_rscu");
	system("rm $out_genes");
	system("perl nfeatures_to_features.pl $tmp_out_features $tmp_out_rscu $final_features");
	system("rm $tmp_out_features");
	system("rm $tmp_out_rscu");
	system("python3 E2_E6_models/E2_E6_load_model.py $final_features $model");
	system("rm $final_features");
	system("rm -r $prefix");
}
elsif($input_type=~m/all_genes/){
	$out_genes="$input_genome";
	@input_genome1=split('\.',$input_genome);
	$prefix=$input_genome1[0];
	$final_features= "$prefix".'_features.csv';
	$tmp_out_features='HPV_features.tsv';
	$tmp_out_rscu='HPV_rscu.tsv';
	#system("perl genome_to_gene.pl $input_genome $out_genes");
	system("perl nfeatures.pl $out_genes $tmp_out_features");
	system("perl codon_analysis.pl $out_genes $tmp_out_rscu");
	system("rm $out_genes");
	system("rm -f $input_genome");
	system("perl nfeatures_to_features.pl $tmp_out_features $tmp_out_rscu $final_features");
	system("rm $tmp_out_features");
	system("rm $tmp_out_rscu");
	system("python3 15f_five_genes/15f_load_model.py $final_features $model");
	system("rm $final_features");
	#system("rm -r $prefix");
}
elsif($input_type=~m/E2_E6/){
	$out_genes="$input_genome";
	@input_genome1=split('\.',$input_genome);
	$prefix=$input_genome1[0];
	$final_features= "$prefix".'_features.csv';
	$tmp_out_features='HPV_features.tsv';
	$tmp_out_rscu='HPV_rscu.tsv';
	#system("perl genome_to_gene.pl $input_genome $out_genes");
	system("perl nfeatures.pl $out_genes $tmp_out_features");
	system("perl codon_analysis.pl $out_genes $tmp_out_rscu");
	system("rm $out_genes");
	system("rm -f $input_genome");
	system("perl nfeatures_to_features.pl $tmp_out_features $tmp_out_rscu $final_features");
	system("rm $tmp_out_features");
	system("rm $tmp_out_rscu");
	system("python3 E2_E6_models/E2_E6_load_model.py $final_features $model");
	system("rm $final_features");
	#system("rm -r $prefix");
}
elsif($input_type=~m/E6/){
	$out_genes="$input_genome";
	@input_genome1=split('\.',$input_genome);
	$prefix=$input_genome1[0];
	$final_features= "$prefix".'_features.csv';
	$tmp_out_features='HPV_features.tsv';
	$tmp_out_rscu='HPV_rscu.tsv';
	#system("perl genome_to_gene.pl $input_genome $out_genes");
	system("perl nfeatures.pl $out_genes $tmp_out_features");
	system("perl codon_analysis.pl $out_genes $tmp_out_rscu");
	system("rm -f $input_genome");
	system("rm -f $out_genes");
	system("perl nfeatures_to_features.pl $tmp_out_features $tmp_out_rscu $final_features");
	system("rm -f $tmp_out_features");
	system("rm -f $tmp_out_rscu");
	system("python3 E6_models/load_E6_model.py $final_features $model");
	system("rm $final_features");
	#system("rm -r $prefix");
}
