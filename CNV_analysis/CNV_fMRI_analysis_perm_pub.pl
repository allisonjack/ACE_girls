#! perl -w

use strict;
use warnings;
use autodie;
use Getopt::Long;
use POSIX;
use Cwd qw(abs_path);
use List::Util qw(shuffle);
use Data::Dumper;


my ($affected, $type, $gender, $sex, $permutations, $brain_region, $time, $exclusion_list, $number_for_groups, $start_from_pulling_list_CNVs_genes, $start_after_CNV_files_made, $start_after_exonarray_filtered, $start_from_new_regions_picked);
GetOptions( 
	'affected' => \$affected,
	'ptype=s' => \$type,
	'gender' => \$gender,
	'sex=s' => \$sex,
	'ipermutations=s' => \$permutations,	
	'brain=s' => \$brain_region,
	'time=s' => \$time,
	'zexclude=s' => \$exclusion_list,
	'numbers=s' => \$number_for_groups,
	'xlist=s' => \$start_from_pulling_list_CNVs_genes,
	'yfiles' => \$start_after_CNV_files_made,
	'wexonarray' => \$start_after_exonarray_filtered,
	'vnewregions' => \$start_from_new_regions_picked);

my $USAGE = "

	#================================================================================================
	### DESCRIPTION: To take an input file with samples and sizes of CNVs and do a permutation test for list of genes of interest
	### and perform permutations to compare and see if significant
	### Please provide the input data file, the means R script and select the options you wish to compare
	### EX: perl CNV_fMRI_analysis_perm_pub.pl <input_data_file><criteria_basename><output_name><exonarray_metadata><columns><genecode> <calculate_means.R><extract_sample_values_info2.R> -p proband -g -t 3,4,5 -b STR_R,STR_L
	### EX2: perl CNV_fMRI_analysis_perm_pub.pl <input_data_file><criteria_basename><output_name><exonarray_metadata><columns><genecode><calculate_means.R><extract_sample_values_info2.R> -p proband,control -t 3,4,5 -b STR_R,STR_L
	### EX3: perl CNV_fMRI_analysis_perm_pub.pl <input_data_file><FPE><output_name><exonarray_metadata><columns><genecode><calculate_means.R><extract_sample_values_info2.R> -p affected -g -y -i 10 -t 3,4,5 -b STR_R,STR_L
	### OUTPUT: Will provide a new text file containing the p-value of the original samples and if permutations their p-values 
	### OPTION: -a to compare data by affected status, -p to give specific phenotypes to compare (ex. proband, sibling, control) 
	### OPTION2: -g to compare data by gender, -s to give a specific gender to look at (ie. female)
	### OPTION3:  -i number of iterations to perform, -y to start after the different CNV files made all with same basename
	### OPTION4: -w to start from after the exonarray data has been filtered for CNV genes
	### OPTION5: -v to start from picking new brain regions of interest and time periods to find genes expressed in
	### OPTION6: -b give brain regions of interest (sep by , :ex. STR_R,DFC_R), -t give time regions of interest (sep by , or - range)
	### OPTION7: -z give list of samples to exclude if have any, 
	### OPTION8: -n (to give number to select for each group sep by , and followed by mean difference total size)
	#==================================================================================================\n";

#END_USAGE

die $USAGE if (@ARGV != 8);

my $input_data = $ARGV[0];  ## file containing header and list of IDs, Gender, Phenotype, CNV Size, Genes
my $criteria_name = $ARGV[1];
my $output_name = $ARGV[2];
my $exonarray_metadata = $ARGV[3]; ##exonarray data file
my $columns_data = $ARGV[4];  ##exonarray map of columns info
my $genecode = $ARGV[5];  ##exonarray map of gene codes
my $rscript_file = $ARGV[6];  ##calculate_means.R
my $rscript_extract_info = $ARGV[7];  ##extract_sample_values_info2.R

if (!-r $input_data) {die "Unable to open input file $input_data!!!\n";}

use Cwd qw(cwd); ##returns current working directory even if a symbolic link
my $dir = cwd;
use Cwd qw(abs_path);

my $logfile = "${output_name}_logfile.txt";
my $LOG1;
if (-e $logfile) {open $LOG1, '>>', $logfile; print $LOG1 "\n\n";}
else {open $LOG1, '>', $logfile;}

print "Currently running program: ".abs_path($0)."\n";
print $LOG1 "Currently running program: ".abs_path($0)."\n";
my $original_starttime = time;
my $starttime = localtime();
print "Starttime is $starttime.\n";
print $LOG1 "Startime is $starttime.\n";
my $abspathin = &abs_path($0);
my $programpath = substr ($abspathin, 0, (rindex($abspathin, '/'))+1);
print "Folder containing programs are $programpath.\n";
my $basename = substr ($input_data, 0, (rindex($input_data, '.')));
print "My input file basename is $basename.\n";

my @exon_regions_avail = ("OFC_R","OFC_L","DFC_R","DFC_L","VFC_R","VFC_L","MFC_R","MFC_L","M1C_R","M1C_L","S1C_R","S1C_L","IPC_R","IPC_L","A1C_R","A1C_L","STC_R","STC_L","ITC_R","ITC_L","V1C_R","V1C_L","HIP_R","HIP_L","AMY_R","AMY_L","STR_R","STR_L","MD_R","MD_L","CBC_R","CBC_L");
my $num_brain_regions_avail = scalar @exon_regions_avail;
if ($brain_region) {
	if ($brain_region =~ /,/){
		my @regions_check = split (",",$brain_region);
		foreach my $check (@regions_check) {
			if (grep (/$check/,@exon_regions_avail)) {next;}
			else {die "At least one brain region you gave is not an option on the exonarray data!!!\nHere are the options: @exon_regions_avail\n\n";}
			}
		}
	else {if (!grep (/$brain_region/,@exon_regions_avail)) {die "The brain region you gave is not an option on the exonarray data!!!\nHere are the options: @exon_regions_avail\n\n";}}
	}
my @times;	
if ($permutations && !$time) {die "You must enter time periods to choose from for the random regions!!!\n";}
if ($time) {
	if ($time =~ /,/) {
		@times = split ',',$time;}
	elsif ($time =~ /-/) {
		my ($first_num,$second_num) = split (/-/,$time);
		my $c=0;
		for (my $n=$first_num; $n<=$second_num; $n++) {
			$times[$c] = ($n);
			$c++;
			}
		}
	else {
		$times[0] = $time;}
	}

if ($type) {print "You want to compare $type.\n\n";}
if ($affected) {print "You want to compare affected status.\n\n";}
if ($gender) {print "You want to compare the genders.\n\n";}
if ($sex) {print "You want to compare $sex.\n\n";}

my $type_column;
my $affected_column;
my $type_avail;
my @types_avail;
my $gender_column;
my $gender_avail;
my @gender_avail;
my $gene_column;
$gene_column = qx {head -1 $input_data | awk '{for (i=1;i<=NF;i++) if(\$i ~ /Gene/){print i;exit}}'}; #which column has gene names
chomp $gene_column;
if (!$gene_column) {die "The file you gave must contain a column for the Genes found in the CNVs!\n";}
my ($id_column) = qx {head -1 $input_data | awk '{for (i=1;i<=NF;i++) if(\$i ~/ID/){print i;exit}}'}; #which column has the IDs
chomp $id_column;
my ($id_col_array) = $id_column -1;
my (@ids, $ids);
$ids = qx{awk -F"\t" 'NR>1 {print \$${id_column}}' $input_data | sort | uniq}; #collect IDs from file
chomp $ids;
@ids = split ("\n", $ids); #array of IDs available
my ($num_ids) = scalar @ids;
my ($num_female, $num_male, $num_affected, $num_unaffected);
my @excluded;
my @final_criteria_genes=();
my @final_criteria_positive_genes=();
my $regions_avail;
my @regions_avail;
my ($num_gp1,$num_gp2,$numbers_mean);
my @number_groups;
if ($number_for_groups) {@number_groups = split ',', $number_for_groups; if (scalar @number_groups != 3 || ($number_groups[0] != int($number_groups[0])) || ($number_groups[1] != int($number_groups[1]))) {die print "Can't have more or less than 2 numbers for comparison (1 for each group) followed by the mean for total CNVs.\n\n";} else {$num_gp1 = $number_groups[0]; $num_gp2 = $number_groups[1]; $numbers_mean = $number_groups[2]}}

my ($input_data_CNV_genes_merge,$input_data_all_CNV_genes,$input_genes_exonarray_data,$input_genes_exonarray_list,$input_CNVs_exonarray,$input_genes_criteria,$input_CNVs_criteria);
my ($input_data_CNVs_merge_ids, $input_data_CNV_exonarray_merge_ids, $input_data_CNV_criteria_merge_ids);

if ($start_after_CNV_files_made) { ##-y
	$input_data_CNV_genes_merge = "${basename}_all_CNV_genelist_merged.txt";#file of CNV genes per individual
	$input_data_all_CNV_genes = "${basename}_CNV_genes_to_get_data_for.txt";#list of the CNV genes to look in exonarray metadata for
	$input_genes_exonarray_data = "${basename}_genes_metadata.txt"; #output name for running pull_exonarray giving genes available and their metadata - will be a file of exonarray data for only those genes found in the CNVs
	$input_genes_exonarray_list = "${basename}_genes_metadata_final_genelist.txt";#output name for pull_exonarray list of genes on exonarray
	$input_CNVs_exonarray = "${basename}_exonarray_data.txt"; #original list selected for only those CNVs with genes on exonarray
	$input_genes_criteria = "${basename}_${criteria_name}_criteria_genes.txt"; #name for pulling out exonarray CNVs within criteria brain regions and times
	$input_CNVs_criteria = "${basename}_${criteria_name}_criteria_exonarray_data.txt"; #exonarray CNVs selected for only the criteria regions of interest
	if (!-e $input_genes_exonarray_data || !-e $input_genes_exonarray_list || !-e $input_CNVs_exonarray || !-e $input_genes_criteria || !-e $input_CNVs_criteria) {die "Unable to find one of the expected files!!!\n";}
	}
elsif ($start_after_exonarray_filtered) { ##-w
	if (!$brain_region || !$time) {die "At least one brain region and one time period needs to be given as an option to continue!!!";}
	$input_data_CNV_genes_merge = "${basename}_all_CNV_genelist_merged.txt";#file of CNV genes per individual
	$input_data_all_CNV_genes = "${basename}_CNV_genes_to_get_data_for.txt";#list of the CNV genes to look in exonarray metadata for
	$input_genes_exonarray_data = "${basename}_genes_metadata.txt"; #output name for running pull_exonarray giving genes available and their metadata
	$input_genes_exonarray_list = "${basename}_genes_metadata_final_genelist.txt";#output name for pull_exonarray list of genes on exonarray
	if (!-e $input_genes_exonarray_data || !-e $input_genes_exonarray_list) {die "Unable to find one of the expected files!!!\n";}
	if (!-r $input_genes_exonarray_data) {die "Unable to read exonarray data file $input_genes_exonarray_data!!!  Check that it exists!!!\n";}
	
	print "Making new file containing only the CNVs per individual that are in the exonarray data.\n";
	$input_CNVs_exonarray = "${basename}"; #original list selected for only those CNVs with genes on exonarray
	&extract_gene_list($input_data, $input_genes_exonarray_list, $input_CNVs_exonarray);
	$input_CNVs_exonarray = $input_CNVs_exonarray."_exonarray_data.txt";
	
	print "Finding which genes involved in CNVs are positively expressed for the regions of interest.\n";
	print "Inputing the options for time $time and brain regions of interest $brain_region into the program.\n\n";
	$input_genes_criteria = "${basename}_${criteria_name}_criteria"; #name for pulling out exonarray CNVs within criteria brain regions and times
	&criteria_expression_median_genes($input_genes_exonarray_data, $input_genes_criteria, $time, $brain_region);
	$input_genes_criteria = $input_genes_criteria."_genes.txt";
	
	print "Making new file containing only the CNVs per individual contain genes that are positive for the regions of interest.\n";
	$input_CNVs_criteria = "${basename}_${criteria_name}_criteria"; #exonarray list selected for only the criteria regions of interest
	&extract_gene_list($input_CNVs_exonarray, $input_genes_criteria, $input_CNVs_criteria);
	$input_CNVs_criteria = $input_CNVs_criteria."_exonarray_data.txt";
	}
elsif ($start_from_new_regions_picked) { ##-v  to use same CNV list and exonarray extracted file and select new regions of interest
	if (!$brain_region || !$time) {die "At least one brain region and one time period needs to be given as an option to continue!!!";}
	$input_data_CNV_genes_merge = "${basename}_all_CNV_genelist_merged.txt";#file of CNV genes per individual
	$input_data_all_CNV_genes = "${basename}_CNV_genes_to_get_data_for.txt";#list of the CNV genes to look in exonarray metadata for
	$input_genes_exonarray_data = "${basename}_genes_metadata.txt"; #output name for running pull_exonarray giving genes available and their metadata
	$input_genes_exonarray_list = "${basename}_genes_metadata_final_genelist.txt";#output name for pull_exonarray list of genes on exonarray
	$input_CNVs_exonarray = "${basename}_exonarray_data.txt"; #original list selected for only those CNVs with genes on exonarray
	if (!-e $input_genes_exonarray_data || !-e $input_genes_exonarray_list || !-e $input_CNVs_exonarray) {die "Unable to find one of the expected files!!!\n";}
	if (!-r $input_genes_exonarray_data) {die "Unable to read exonarray data file $input_genes_exonarray_data!!!  Check that it exists!!!\n";}
	
	print "Finding which genes involved in CNVs are positively expressed for the regions of interest.\n";
	print "Inputing the options for time $time and brain regions of interest $brain_region into the program.\n\n";
	$input_genes_criteria = "${basename}_${criteria_name}_criteria"; #name for pulling out exonarray CNVs within criteria brain regions and times
	&criteria_expression_median_genes($input_genes_exonarray_data, $input_genes_criteria, $brain_region, $time);
	$input_genes_criteria = $input_genes_criteria."_genes.txt";
	
	print "Making new file containing only the CNVs per individual contain genes that are positive for the regions of interest.\n";
	$input_CNVs_criteria = "${basename}_${criteria_name}_criteria"; #exonarray list selected for only the criteria regions of interest
	&extract_gene_list($input_CNVs_exonarray,$input_genes_criteria,$input_CNVs_criteria);
	$input_CNVs_criteria = $input_CNVs_criteria."_exonarray_data.txt";
	}
else {
	if (!$brain_region || !$time) {die "At least one brain region and one time period needs to be given as an option to continue!!!";}
	print "Making merging file of all genes in CNVs per individual.\n";
	my $temp_CNV_to_merge = "${basename}_prog_temp_to_merge.txt"; #temp file for merging CNV genes
	qx {awk -F"\t" 'BEGIN{OFS=FS} {print \$${id_column},\$${gene_column}}' $input_data | sed 's/\"//g' | sed 's/\r//g' > $temp_CNV_to_merge};
	my ($num_lines_to_merge) = qx {wc -l < $temp_CNV_to_merge};
	chomp $num_lines_to_merge;
	print "There are $num_lines_to_merge number of CNVs in the original files.\n";
	print $LOG1 "There are $num_lines_to_merge number of CNVs in the original files.\n";
	$input_data_CNV_genes_merge = "${basename}_all_CNV_genelist_merged.txt"; #file of CNV genes per individual
	qx {head -1 $temp_CNV_to_merge > $input_data_CNV_genes_merge};
	qx {awk -F"\t" 'NR>1 {a[\$1]=(\$1 in a? a[\$1]",":"")\$2}END{for(i in a)print i,a[i]}' OFS="\t" $temp_CNV_to_merge | sort >> $input_data_CNV_genes_merge}; #put into 1 line per individual for all genes in CNVs
	qx {rm $temp_CNV_to_merge}; #remove temp file for merging genes
	my ($num_lines_ind_orig) = qx {wc -l < $input_data_CNV_genes_merge};
	chomp $num_lines_ind_orig;
	print "There are $num_lines_ind_orig number of individuals with CNVs merged from the original file.\n";
	print $LOG1 "There are $num_lines_ind_orig number of individuals with CNVs merged from the original file.\n";
	print "Making file for list of genes from CNVs to look in exonarray data for.\n";
	my $temp_all_CNV_genes = "${basename}_temp_all_CNV_genes.txt";
	$input_data_all_CNV_genes = "${basename}_CNV_genes_to_get_data_for.txt"; #list of the CNV genes to look in exonarray metadata for
	qx {awk -F"\t" 'NR>1 {print \$2}' $input_data_CNV_genes_merge > $temp_all_CNV_genes};
	qx {sed 's/,/\\n/g' $temp_all_CNV_genes | sort | uniq > $input_data_all_CNV_genes};
	my ($num_list_all_CNV_genes) = qx {wc -l < $input_data_all_CNV_genes};
	chomp $num_list_all_CNV_genes;
	print "There are $num_list_all_CNV_genes number of CNV genes from the original file.\n";
	print $LOG1 "There are $num_list_all_CNV_genes number of CNV genes from the original file.\n";
	$input_genes_exonarray_data = "${basename}_genes_metadata"; #output name for running pull_exonarray giving genes available and their metadata
	$input_genes_exonarray_list = "${basename}_genes_metadata_final_genelist.txt"; #output name for pull_exonarray list of genes on exonarray

	print "Pulling out the exonarray data for genes involved in CNVs of interest.\n";
	&pulling_out_exonarray_data($input_data_all_CNV_genes, $exonarray_metadata, $columns_data, $genecode, $input_genes_exonarray_data);
	#now have file of gene exonarray data that were involved in the CNVs and the list of which genes they are
	$input_genes_exonarray_data = $input_genes_exonarray_data.".txt";

	print "Making new file containing only the CNVs per individual that are in the exonarray data.\n";
	$input_CNVs_exonarray = "${basename}"; #original list selected for only those CNVs with genes on exonarray
	&extract_gene_list($input_data, $input_genes_exonarray_list, $input_CNVs_exonarray);
	$input_CNVs_exonarray = $input_CNVs_exonarray."_exonarray_data.txt";

	print "Finding which genes involved in CNVs are positively expressed for the regions of interest.\n";
	print "Inputing the options for time $time and brain regions of interest $brain_region into the program.\n\n";
	$input_genes_criteria = "${basename}_${criteria_name}_criteria"; #name for pulling out exonarray CNVs within criteria brain regions and times
	&criteria_expression_median_genes($input_genes_exonarray_data, $input_genes_criteria, $brain_region, $time);
	$input_genes_criteria = $input_genes_criteria."_genes.txt";
	
	print "Making new file containing only the CNVs per individual contain genes that are positive for the regions of interest.\n";
	$input_CNVs_criteria = "${basename}_${criteria_name}_criteria"; #exonarray list selected for only the criteria regions of interest
	&extract_gene_list($input_CNVs_exonarray,$input_genes_criteria,$input_CNVs_criteria);
	$input_CNVs_criteria = $input_CNVs_criteria."_exonarray_data.txt";
}

print "Original CNV file is $input_data.\n";
print "The exonarray data for available genes found in the CNVs is in $input_genes_exonarray_data.\n";
print "The list of CNV genes that were available in the exonarray is in $input_genes_exonarray_list.\n";
print "The list of CNVs containing genes from the exonarray are in $input_CNVs_exonarray.\n";
print "The list of CNV genes that were expressed in the brain regions of interest are in $input_genes_criteria.\n";
print "The list of CNVs containing genes in the brain regions of interest are in $input_CNVs_criteria.\n\n";

print $LOG1 "Original CNV file is $input_data.\n";
print $LOG1 "The exonarray data for available genes found in the CNVs is in $input_genes_exonarray_data.\n";
print $LOG1 "The list of CNV genes that were available in the exonarray is in $input_genes_exonarray_list.\n";
print $LOG1 "The list of CNVs containing genes from the exonarray are in $input_CNVs_exonarray.\n";
print $LOG1 "The list of CNV genes that were expressed in the brain regions of interest are in $input_genes_criteria.\n";
print $LOG1 "The list of CNVs containing genes in the brain regions of interest are in $input_CNVs_criteria.\n";
$starttime = localtime();
print $LOG1 "$starttime.\n";

my $positive_genes = qx {awk '{print \$0}' $input_genes_criteria}; 
chomp $positive_genes; 
my @criteria_positive_genes = split("\n",$positive_genes); 
if ($criteria_positive_genes[0] =~ /Gene/i) {shift @criteria_positive_genes;}
foreach my $pos (@criteria_positive_genes) {
	$pos =~ s/\"//g;}
print "There are ".scalar @criteria_positive_genes." positive genes from the CNVs in the criteria regions.\n";
print $LOG1 "There are ".scalar @criteria_positive_genes." positive genes from the CNVs in the criteria regions.\n";

my ($type_col_array, $aff_col_array, $gender_col_array);
my $cutoff = 6; #cutoff for positive expression in exonarray data
my %samples;
my %samples_exon;
my %samples_criteria;
my (%samples_total,%samples_exon_total,%samples_crit_total,%samples_critlist_total);
my %exonarray_data;

print "Finding columns of interest in $input_data.\n\n";
	$type_column = qx {head -1 $input_data | awk '{for (i=1;i<=NF;i++) if(\$i ~ /Relation|Phenotype|pheno/){print i}}'};
	chomp $type_column;
	print "Type column is column number $type_column.\n";
	$type_avail = qx{awk -F"\t" 'NR>1 {print \$${type_column}}' $input_data | sort | uniq}; 
	chomp $type_avail; 
	@types_avail = split ("\n",$type_avail);
	$affected_column = qx {head -1 $input_data | awk -F"\t" '{for(i=1;i<=NF;i++)if(\$i~/Affected|Unaffected|Type|status/){print i;exit}}'};
	chomp $affected_column;
	print "The Affected status column is $affected_column.\n";
	$gender_column = qx {head -1 $input_data | awk '{for (i=1;i<=NF;i++) if(\$i ~ /Sex|Gender/){print i;exit}}'};
	chomp $gender_column;
	print "My genders are in column $gender_column.\n";

if ($exclusion_list) {
	print "You have provided a list of samples to be excluded.\n\n";
	my $exID_col = qx {head -1 $exclusion_list | awk '{for (i=1;i<=NF;i++) if(\$i ~ /ID|id/){print i}}'};
	chomp $exID_col;
	my $exclude = qx {awk -F"\t" 'NR>1 {print \$${exID_col}}' $exclusion_list | sort | uniq};
	chomp $exclude;
	@excluded = split ("\n",$exclude);
	@excluded = grep(s/\s*$//g, @excluded);
	print "The samples that need to be excluded are @excluded.\n\n";
	print $LOG1 "You included a list of samples to exclude.  The samples that need to be excluded are @excluded.\n\n";
	}

$type_col_array = $type_column -1;	
my @gender_headers = ('Gender','gender','Sex','sex');

#open input data file to look at each row
open INFILE , '<' , $input_data or die "Could not open file $input_data\n$!";
my $headers = qx {head -1 $input_data};
chomp $headers;
$headers =~ s/[\n\r]+//g;
my @headers = split("\t",$headers); 
my @output_headers;
@output_headers = @headers;
if ($affected_column eq "") {splice @output_headers, ($type_col_array+1), 0, 'Sample_Type';$aff_col_array = $type_col_array +1;$gender_col_array = $gender_column;}
else {$aff_col_array = $affected_column -1;$gender_col_array = $gender_column-1;}
my $headers_count = scalar @output_headers;
print "@output_headers.\n";
my ($size_col) = grep {$headers[$_] =~ /size|Size/} 0..$#headers; 
my ($size_col_array) = grep {$output_headers[$_] =~ /size|Size/} 0..$#output_headers; 
my ($gene_in_col) = grep {$headers[$_] =~ /Gene|gene/} 0..$#headers;
my ($gene_col_array) = grep {$output_headers[$_] =~ /Gene|gene/} 0..$#output_headers;
print "My CNV sizes are in column $size_col.\n";
print "My CNV sizes is in header position $size_col_array.\n";
#$gender_col_array = grep {$output_headers[$_] =~ /Gender|gender|Sex|sex/} 0..$#output_headers;
#$gender_col_array = grep {$output_headers[$_] =~ @gender_headers} 0..$#output_headers;
print "My genders are in header position $gender_col_array.\n";
print "My genes are in header position $gene_col_array.\n";
print "The header for genes is $output_headers[$gene_col_array].\n";
	
print "\nSeeing what terminology is used for the genders.\n\n";
my @male_genders = ('M','m','Male','male','MALE');
my ($male_gender,$female_gender); #finding which annotation used in file for genders and capturing them in these variables
my @female_genders= ('F','f','Female','female','FEMALE');
$gender_avail = qx{awk -F"\t" 'NR>1 {print \$${gender_column}}' $input_data | sort | uniq}; 
chomp $gender_avail; 
@gender_avail = split ("\n",$gender_avail);
foreach my $genders (@gender_avail){
	foreach my $males (@male_genders){
		if ($genders eq $males) {$male_gender = $genders;print "The male gender used in the file is $male_gender.\n";}}
	foreach my $females (@female_genders){
		if ($genders eq $females) {$female_gender = $genders;print "The female gender used in the file is $female_gender.\n";}}
		}
print "\n\n";
my @genders;
my $gender_line;

my @affected = ("proband","Proband","PROBAND","Affected Sibling","affected sibling","Affected sibling","AFFECTED SIBLING");
my @unaffected = ("Sibling","sibling","SIBLING","Unaffected Sibling","Unaffected sibling","unaffected sibling","UNAFFECTED SIBLING","control","Control","CONTROL","Control sibling","Control Sibling","control sibling","CONTROL SIBLING");
my @affect_status_options = ("Affected","AFFECTED","affected","affect","Unaffected","UNAFFECTED","unaffected","unaffect");
my @a_status = ("Affected","AFFECTED","affected","affect","A","ASD");
my @u_status = ("Unaffected","UNAFFECTED","unaffected","unaffect","U","TD","US");
my ($affecteds, $unaffected, $proband, $asib, $usib, $control);
my $pheno_avail = qx{awk -F"\t" 'NR>1 {print \$${type_column}}' $input_data | sort | uniq};
chomp $pheno_avail;
my @pheno_avail = split ("\n",$pheno_avail);
foreach my $pheno (@pheno_avail) {
	if ($pheno =~ /pro/i) {$proband = $pheno;}
	elsif ($pheno =~/con/i && $pheno !~ /sib/i) {$control = $pheno;}
	elsif ($pheno =~ /^aff/i) {$asib = $pheno;}
	elsif ($pheno =~/^un/i || $pheno =~ /^sib/i) {$usib = $pheno;}
	}
my ($aff_avail, @aff_avail);
if ($affected_column ne "") {
	$aff_avail = qx{awk -F"\t" 'NR>1 {print \$${affected_column}}' $input_data | sort | uniq};
	chomp $aff_avail; 
	@aff_avail = split ("\n",$aff_avail);
	foreach my $aff (@aff_avail) {
		if (grep (/$aff/, @u_status)) {$unaffected = $aff;}
		else {$affecteds = $aff;}
			}
		}
else {$unaffected = "Unaffected";$affecteds = "Affected";}
#print "These are found in the file: $proband, $asib, $usib, $control, $affecteds, $unaffected.\n";
my @types;
my $type_count;
my $type_col1;
my $type_col2;
my $type1;
my $type2;
my $num_type1;
my $num_type2;


#determine what option to run
die "No option given to compare\n" unless ($affected || $gender || $type || $sex);
if ($affected) {
	#@types = ('Affected','Unaffected');
	push @types, $affecteds;
	push @types, $unaffected;
	print $LOG1 "Looking at affected status for comparison.\n";
	print "Type1 is $types[0] and Type2 is $types[1].\n";
	$type1 = $types[0];
	$type2 = $types[1];
	}
if ($type) {
	print "You chose to compare type/types.\n\n";
	@types = split(',',$type); 
	print "To compare $type.\n";
	$type_count = scalar @types;
	if ($type_count == 1) {
		if (grep {$types[0] eq $_} @affected) {
			if ($types[0] =~ /pro/i) {$type1 = $proband;}
			else {$type1 = $asib; }
			$type_col1 = $type_col_array;
			}
		elsif (grep {$types[0] eq $_} @unaffected) {
			if ($types[0] =~ /^con/i) {$type1 = $control;}
			else {$type1 = $usib; }
			$type_col1 = $type_col_array;
			}
		elsif (grep {$types[0] eq $_} @affect_status_options) {
			if ($types[0] =~ /^un/i) {$type1 = $unaffected;} 
			else {$type1 = $affecteds;}
			$type_col1 = $aff_col_array;
			}
		print "Looking at type $type1 in position $type_col1 of the headers array.\n";
		print $LOG1 "Looking at type $type1 in position $type_col1 of the headers array.\n";
		}
	elsif ($type_count == 2) {
		for (my $c=1; $c<=2; $c++) {
			if ($c == 1) {
				if (grep {$types[$c] eq $_} @affected) {
					if ($types[$c] =~ /pro/i) {$type1 = $proband;}
					else {$type1 = $asib; }
					$type_col1 = $type_col_array;
					}
				elsif (grep {$types[$c] eq $_} @unaffected) {
					if ($types[$c] =~ /^con/i) {$type1 = $control;}
					else {$type1 = $usib;}
					$type_col1 = $type_col_array;
					}
				elsif (grep {$types[$c] eq $_} @affect_status_options) {
					if ($types[$c] =~ /^un/i) {$type1 = $unaffected;}
					else {$type1 = $affecteds;}
					$type_col1 = $aff_col_array;
					}
				}
			else {
				if (grep {$types[$c] eq $_} @affected) {
					if ($types[$c] =~ /pro/i) {$type2 = $proband;}
					else {$type2 = $asib;}
					$type_col2 = $type_col_array;
					}
				elsif (grep {$types[$c] eq $_} @unaffected) {
					if ($types[$c] =~ /^con/i) {$type2 = $control;}
					else {$type2 = $usib;}
					$type_col2 = $type_col_array;
					}	
				elsif (grep {$types[$c] eq $_} @affect_status_options) {
					if ($types[$c] =~ /^un/i) {$type2 = $unaffected;}
					else {$type2 = $affecteds;}
					$type_col2 = $aff_col_array;
					}
				}
			}
		print "Looking at type $type1 and $type2 in position $type_col1 and position $type_col2 of the headers array.\n";
		print $LOG1 "Looking at type $type1 and $type2 in position $type_col1 and position $type_col2 of the headers array.\n";
		}
	else {die print "You cannot choose more than 2 types to compare.\n\n";}
	if ($number_for_groups) {
		my @temp_num = split(',',$number_for_groups);
		$num_type1 = $temp_num[0]; 
		$num_type2 = $temp_num[1];
		print "You want to look at $num_type1 for type1 and $num_type2 for type2 in the comparisons for numbers.\n";
		print $LOG1 "You want to look at $num_type1 for type1 and $num_type2 for type2 in the comparisons for numbers.\n";
		}
	}

if ($gender) {
	@genders = ($male_gender,$female_gender); 
	print "You want to compare the genders!\n\n";
	print $LOG1 "You want to compare the genders!\n\n";
	}
if ($sex) {
	my @sexes = split(',',$sex); 
	my ($num_sex) = scalar @sexes; 
	if ($number_for_groups) {
		foreach my $sex (@sexes){
			if (grep {$sex eq $_} @gender_avail) {push @genders, $sex;}
			elsif (grep {$sex eq $_} @male_genders) {push @genders, $male_gender;}
			elsif (grep {$sex eq $_} @female_genders) {push @genders, $female_gender;}
			else {die "This $sex is not available on the input file.\n";} 
			}
			my @temp_num = split(',',$number_for_groups);
			if (grep (/$genders[0]/, @male_genders)) {
				$num_male = $temp_num[0]; 
				$num_female = $temp_num[1]; 
				$num_gp1 = $temp_num[1]; 
				$num_gp2 = $temp_num[0];
				}
			else {
				my @temp_num = split(',',$number_for_groups); 
				$num_male = $temp_num[1]; 
				$num_female = $temp_num[0]; 
				$num_gp1 = $temp_num[0]; 
				$num_gp2 = $temp_num[1];
				}
		$gender = 1;
		}
	elsif ($num_sex != 1) {print "Please select gender option to look at both sexes.\n\n"; die;} 
	else {
		if (grep {$sex eq $_} @gender_avail) {push @genders, $sex;}
		elsif (grep {$sex eq $_} @male_genders) {push @genders, $male_gender;}
		elsif (grep {$sex eq $_} @female_genders) {push @genders, $female_gender;}
		else {die "This $sex is not available on the input file.\n";}
		print "To only look at $sex.\n";	
		print $LOG1 "You want to only look at $sex.\n\n";	
		}
	}
	
if ($type && $gender) {print "You wish to compare the genders of type $type.\n\n";}
elsif ($sex && $affected) {print "You wish to compare the affected status for gender $sex.\n\n";}
elsif ($sex && $type) {print "You wish to compare the types $type for the gender $sex.\n\n";}
elsif ($affected) {print "You wish to compare the affected status for the cohort.\n\n";}
elsif ($type) {print "You wish to compare the types $type for the cohort.\n\n";}


#looking at the different CNV files
my $total_lines_orig = qx {wc -l < $input_data};
chomp $total_lines_orig;
my $total_lines_exon = qx {wc -l < $input_CNVs_exonarray};
chomp $total_lines_exon;
my $total_lines_criteria = qx {wc -l < $input_CNVs_criteria};
chomp $total_lines_criteria;
my ($count_line1, $count_line2, $count_line3) = 2;

$starttime = localtime();
print "Making hashes to contain the different CNV datas.\n\n";
print $LOG1 "$starttime:   Making hashes to contain the different CNV datas.\n";
print $LOG1 "$starttime:   Making hash to contain the original CNV data.\n";
#make hash of original CNV file
print "\nMaking hash for the original CNV data.\n";
for (my $i=2;  $i<=$total_lines_orig; $i++){
	if (int($i/50) == ($i/50)) {print "Looking at CNV $i.\n";}
	if ($i==$total_lines_orig) {print "Looking at CNV $i.\n";}
	my $line = qx {sed -n '${i}p' $input_data};
	$line =~ s/[\n\r]+//g;
	#if (int($i/10) == ($i/10)) {print "$line.\n";}
	my @tab = split /\t/,$line;
	for (my $j=0; $j<$headers_count; $j++){
		if ($output_headers[$j] eq "Sample_Type" && ($affected_column eq "" || !$affected_column)) {
			if (grep {$tab[$type_col_array] eq $_} @affected) {$samples{$i}{$output_headers[$j]} = "Affected";}
			else {$samples{$i}{$output_headers[$j]} = "Unaffected";}
			}
		elsif ($affected_column eq "" || !$affected_column) {
			if ($j == $gene_col_array) {$tab[$gene_in_col] =~ s/\"//g; $samples{$i}{$output_headers[$j]} = $tab[$gene_in_col];
				if (int($i/50) == ($i/50)) {print "CNV $i for $output_headers[$j] is *$tab[$gene_in_col]*.\n";}
			}
			elsif ($j < $aff_col_array) {$samples{$i}{$output_headers[$j]} = $tab[$j];
				if (int($i/50) == ($i/50)) {print "CNV $i for $output_headers[$j] is *$tab[$j]*.\n";}
				}
			else {$samples{$i}{$output_headers[$j]} = $tab[$j-1];
				if (int($i/50) == ($i/50)) {print "CNV $i for $output_headers[$j] is *$tab[($j-1)]*.\n";}
				}
			}
		else {
			if ($j == $gene_col_array) {$tab[$j] =~ s/\"//g; $samples{$i}{$output_headers[$j]} = $tab[$j];
				if (int($i/50) == ($i/50)) {print "CNV $i for $output_headers[$j] is *$tab[$j]*.\n";}
			}
			else {$samples{$i}{$output_headers[$j]} = $tab[$j];
				if (int($i/50) == ($i/50)) {print "CNV $i for $output_headers[$j] is *$tab[$j]*.\n"}
				}
			}
		#if (int($i/10) == ($i/10)) {
		#	print "j is **$j**:  output headers is **$output_headers[$j]**:  tab is **$tab[$j]**\n";}
		
		}
	}

$starttime = localtime();
#make hash of exonarray CNV file
print "\nMaking hash for the exonarray available CNV data.\n";
print $LOG1 "$starttime:   Making hash to contain the exonarray available CNV data.\n";
for (my $i=2;  $i<=$total_lines_exon; $i++){
	if (int($i/50) == ($i/50)) {print "Looking at exonarray CNV $i.\n";}
	if ($i==$total_lines_exon) {print "Looking at exonarray CNV $i.\n";}
	my $line = qx {sed -n '${i}p' $input_CNVs_exonarray};
	$line =~ s/[\n\r]+//g;
	my @tab = split /\t/,$line;
	for (my $j=0; $j<$headers_count; $j++){
		if ($output_headers[$j] eq "Sample_Type" && ($affected_column eq "" || !$affected_column)) {
			if (grep {$tab[$type_col_array] eq $_} @affected) {$samples_exon{$i}{$output_headers[$j]} = "Affected";}
			else {$samples_exon{$i}{$output_headers[$j]} = "Unaffected";}
			}
		elsif ($affected_column eq "" || !$affected_column) {
			if ($j == $gene_col_array) {$tab[$gene_in_col] =~ s/\"//g; $samples_exon{$i}{$output_headers[$j]} = $tab[$gene_in_col];
				if (int($i/50) == ($i/50)) {print "Exonarray CNV $i for $output_headers[$j] is *$tab[$gene_in_col]*.\n";}
			}
			elsif ($j < $aff_col_array) {$samples_exon{$i}{$output_headers[$j]} = $tab[$j];
				if (int($i/50) == ($i/50)) {print "Exonarray CNV $i for $output_headers[$j] is *$tab[$j]*.\n";}
				}
			else {$samples_exon{$i}{$output_headers[$j]} = $tab[$j-1];
				if (int($i/50) == ($i/50)) {print "Exonarray CNV $i for $output_headers[$j] is *$tab[$j-1]*.\n";}
				}
			}
		else {
			if ($j == $gene_col_array) {$tab[$j] =~ s/\"//g; $samples_exon{$i}{$output_headers[$j]} = $tab[$j];
				if (int($i/50) == ($i/50)) {print "Exonarray CNV $i for $output_headers[$j] is *$tab[$j]*.\n";}
			}
			else {
				$samples_exon{$i}{$output_headers[$j]} = $tab[$j];
				if (int($i/50) == ($i/50)) {print "Exonarray CNV $i for $output_headers[$j] is *$tab[$j]*.\n";}
				}
			}
		}
	}	

$starttime = localtime();
#make hash of regions of interest CNV file
print "\nMaking hash for the regions of interest available CNV data.\n\n";
print $LOG1 "$starttime:   Making hash to contain the regions of interest available CNV data.\n";
for (my $i=2;  $i<=$total_lines_criteria; $i++){
	if (int($i/50) == ($i/50)) {print "Looking at regions CNV $i.\n";}
	if ($i==$total_lines_criteria) {print "Looking at regions CNV $i.\n";}
	my $line = qx {sed -n '${i}p' $input_CNVs_criteria};
	$line =~ s/[\n\r]+//g;
	my @tab = split /\t/,$line;
	for (my $j=0; $j<$headers_count; $j++){
		if ($output_headers[$j] eq "Sample_Type" && ($affected_column eq "" || !$affected_column)) {
			if (grep {$tab[$type_col_array] eq $_} @affected) {$samples_criteria{$i}{$output_headers[$j]} = "Affected";}
			else {$samples_criteria{$i}{$output_headers[$j]} = "Unaffected";}
			}
		elsif ($affected_column eq "" || !$affected_column) {
			if ($j == $gene_col_array) {$tab[$gene_in_col] =~ s/\"//g; $samples_criteria{$i}{$output_headers[$j]} = $tab[$gene_in_col];
				if (int($i/50) == ($i/50)) {print "Criteria CNV $i for $output_headers[$j] is *$tab[$gene_in_col]*.\n";}
			}
			elsif ($j < $aff_col_array) {$samples_criteria{$i}{$output_headers[$j]} = $tab[$j];
				if (int($i/50) == ($i/50)) {print "Criteria CNV $i for $output_headers[$j] is *$tab[$j]*.\n";}
				}
			else {$samples_criteria{$i}{$output_headers[$j]} = $tab[$j-1];
				if (int($i/50) == ($i/50)) {print "Criteria CNV $i for $output_headers[$j] is *$tab[$j-1]*.\n";}
				}
			}
		else {
			if ($j == $gene_col_array) {$tab[$j] =~ s/\"//g; $samples_criteria{$i}{$output_headers[$j]} = $tab[$j];
				if (int($i/50) == ($i/50)) {print "Criteria CNV $i for $output_headers[$j] is *$tab[$j]*.\n";}
				}
			else {
				$samples_criteria{$i}{$output_headers[$j]} = $tab[$j];
				if (int($i/50) == ($i/50)) {print "Criteria CNV $i for $output_headers[$j] is *$tab[$j]*.\n";}
				}
			}
		}
	}	
	
if ($affected_column eq "") {$affected_column = $type_column + 1;}	
print "\n\n";

### Names of output files...
$starttime = localtime();
print $LOG1 "$starttime:  Making the different output files.\n";
my $outfile_all= "${output_name}_all_cnvs_permutations.txt"; #picking from all available
print "Making output file $outfile_all.\n";
open my $OUT, '>', $outfile_all;
my $outfile_exon = "${output_name}_exonarray_permutations.txt"; #picking from any on the exonarray 
print "Making output file based on exonarray CNVs only $outfile_exon.\n";
open my $OUTE, '>', $outfile_exon;
my $outCriteria = "${output_name}_criteria_shuffling_permutations.txt"; #label shuffling of criteria groups
print "Making output file based on criteria's permutation shuffling method $outCriteria.\n";
open my $OUTC, '>', $outCriteria;
my $outfile_group = "${output_name}_group_permutations.txt";  #only selected from group of interest ie male/female probands
print "Making output file based on group of interest $outfile_group.\n";
open my $OUTG, '>', $outfile_group;
my $outfile_exon_group = "${output_name}_exon_group_permutations.txt"; #picking from exonarray in group of choice ie male/female probands
print "Making output file based on group of interest on the exonarray $outfile_exon_group.\n";
open my $OUTEG, '>', $outfile_exon_group;

my ($OUTN,$outfile_numbers);
if ($number_for_groups){
	$outfile_numbers = "${output_name}_numbers_permutations.txt"; #picking from any in the regions of interest for specific numbers ie certain number per sex 
	print "Making output file based on numbers of CNVs to choose $outfile_numbers.\n";
	open $OUTN, '>', $outfile_numbers;
	}

my $outfile_final = "${output_name}_permutations_final_results.txt";
print "Making output file for final results from mean difference permutations $outfile_final.\n";
open my $OUTF, '>', $outfile_final;
my $outfile_final_total = "${output_name}_CNVs_total_size.txt";
open my $OUTFT, '>', $outfile_final_total;
my $outfile_final_info = "${output_name}_CNVs_totals_info.txt";
open my $OUTI, '>>', $outfile_final_info;

	
###Start making files
print $OUT "Run\tData file\tComparison\tMean Difference Total Size\n";
print $OUTC "Run\tData file\tComparison\tMean Difference Total Size\n";
print $OUTE "Run\tData file\tComparison\tMean Difference Total Size\n";
print $OUTEG "Run\tData file\tComparison\tMean Difference Total Size\n";
if ($number_for_groups) {
	print $OUTN "Run\tData file\tComparison\tMean Difference Total Size\n";
	}
print $OUTG "Run\tData file\tComparison\tMean Difference Total Size\n";
if ($gender) {
	print $OUTF "File used\tComparison\tMean Difference Total Size\tNum CNVs Female\tNum CNVs Male\tNum Samples Female Total CNVs\tNum Samples Male Total CNVs\n";
	}
else {
	print $OUTF "File used\tComparison\tMean Difference Total Size\tNum CNVs Group1\tNum CNVs Group2\tNum Samples $type1 Total CNVs\tNum Samples $type2 Total CNVs\n";
	}
	my $comparison_text;

my $count =1;

my $temp_outfile = "${output_name}_mean_diff_all.txt";
my $total_outfile = "${output_name}_mean_diff_all_total.txt";
my $temp_outfile1 = "${output_name}_mean_diff_exon.txt";
my $total_exon_outfile = "${output_name}_mean_diff_exon_total.txt";
my $temp_outfile2 = "${output_name}_mean_diff_criteria.txt";
my $total_crit_outfile = "${output_name}_mean_diff_criteria_total.txt";
open my $OUT_temp, '>', $temp_outfile;
open my $OUT_temp1, '>', $temp_outfile1;
open my $OUT_temp2, '>', $temp_outfile2;
open my $OUT_total, '>', $total_outfile;
open my $OUT_total1, '>', $total_exon_outfile;
open my $OUT_total2, '>', $total_crit_outfile;

my $count1;
my $count2;

my (@sample_keys) = keys %samples;
my (@exon_keys) = keys %samples_exon;
my (@criteria_keys) = keys %samples_criteria;
my (@group1,@exon1,@criteria1,@list_gp1);
my ($count_group1,$count_exon1,$count_crit1,$count_list1)=(0,0,0,0);
my (@group2,@exon2,@criteria2,@list_gp2);
my ($count_group2,$count_exon2,$count_crit2,$count_list2)=(0,0,0,0);
my (@criteria_group,@all_group_list);
my ($criteria_count,$all_count_list)=(0,0);
my (@ids_total,@ids_exon_total,@ids_crit_total);
my ($total_crit_mean_diff,$total_mean_diff,$total_exon_mean_diff);
my ($numcnvs1,$numcnvs2)=(0)x2;
my ($totalind1,$totalind2,$totalexonind1,$totalexonind2,$totalcritind1,$totalcritind2)=(0)x6;
my ($info_data_all,$info_data_exon,$info_data_criteria);
my (@check_crit_cnv_id,@check_critlist_cnv_id);

#print "Before going into the selections to extract data, type1 is $type1 and type2 is $type2.\n\n";

#if ($gender || $sex && $affected || $type) {
	if ($gender && $type) {
		print "Comparing the genders for $type.\n";print "Will look for $type1 in $type_col1.\n";print "Will look for gender in $gender_col_array.\n";print "Will extract the CNV size from $size_col_array.\n";
		print $LOG1 "Comparing the genders for $type.\nWill look for $type1 in $type_col1.\nWill look for gender in $gender_col_array.\nWill extract the CNV size from $size_col_array.\n\nLooking at all CNVs.\n";
		print $OUT_temp "Gender\tSize\n";
		print $OUT_temp1 "Gender\tSize\n";
		print $OUT_temp2 "Gender\tSize\n";
		print $OUT_total "Gender\tSize\n";
		print $OUT_total1 "Gender\tSize\n";
		print $OUT_total2 "Gender\tSize\n";
		print "\nLooking at all CNVs.\n";
		foreach my $key (@sample_keys){
			#print "Looking at key $key which is for sample $samples{$key}{$output_headers[$id_col_array]}.\n";
			if (grep {$samples{$key}{$output_headers[$id_col_array]} eq $_} @excluded) {
				#print "The sample $samples{$key}{$output_headers[$id_col_array]} is excluded.\n"; 
				next;} 
			#print "Sample $samples{$key}{$output_headers[$id_col_array]} is not excluded.\n";
			#print "Sample $samples{$key}{$output_headers[$id_col_array]} is type $samples{$key}{$output_headers[$type_col1]}.\n"; 
			if ($samples{$key}{$output_headers[$type_col1]} eq $type1) {
				#print "Sample is $type1.\n";
				if ($samples{$key}{$output_headers[$gender_col_array]} eq $male_gender) {
					push @group2, $key;
					$count_group2+=1;
					#print $OUT_temp "Male\t".$samples{$key}{$output_headers[$size_col_array]}."\n";
					#print "Sample is a male.  CNV size is $samples{$key}{$output_headers[$size_col_array]}\n";
					if (grep {$samples{$key}{$output_headers[$id_col_array]} eq $_} @ids_total) {
						#print "$samples{$key}{$output_headers[$id_col_array]}\n";	
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Size} = $samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Size}+$samples{$key}{$output_headers[$size_col_array]};
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Genes} = $samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Genes}.",".$samples{$key}{$output_headers[$gene_col_array]};
						}
					else {
						#print "$samples{$key}{$output_headers[$id_col_array]}\t";
						push @ids_total, $samples{$key}{$output_headers[$id_col_array]};
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{ID} = $samples{$key}{$output_headers[$id_col_array]};
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Pheno} = $samples{$key}{$output_headers[$aff_col_array]};
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Type} = $samples{$key}{$output_headers[$type_col_array]};
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Gender} = $samples{$key}{$output_headers[$gender_col_array]};
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Size} = $samples{$key}{$output_headers[$size_col_array]};
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Genes} = $samples{$key}{$output_headers[$gene_col_array]};
						$totalind2+=1;
						#print "$totalind2\n";
						}
					}
				if ($samples{$key}{$output_headers[$gender_col_array]} eq $female_gender) {
					push @group1, $key; 
					$count_group1+=1;
					#print $OUT_temp "Female\t".$samples{$key}{$output_headers[$size_col_array]}."\n";
					#print "Sample is a female.\n";
					if (grep {$samples{$key}{$output_headers[$id_col_array]} eq $_} @ids_total) {
						#print "$samples{$key}{$output_headers[$id_col_array]}\n";
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Size} = $samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Size}+$samples{$key}{$output_headers[$size_col_array]};
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Genes} = $samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Genes}.",".$samples{$key}{$output_headers[$gene_col_array]};
						}
					else {
						#print "$samples{$key}{$output_headers[$id_col_array]}\t";
						push @ids_total, $samples{$key}{$output_headers[$id_col_array]};
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{ID} = $samples{$key}{$output_headers[$id_col_array]};
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Pheno} = $samples{$key}{$output_headers[$aff_col_array]};
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Type} = $samples{$key}{$output_headers[$type_col_array]};
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Gender} = $samples{$key}{$output_headers[$gender_col_array]};
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Size} = $samples{$key}{$output_headers[$size_col_array]};
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Genes} = $samples{$key}{$output_headers[$gene_col_array]};
						$totalind1+=1;
						#print "$totalind1\n";
						}
					}
				}
			}
		foreach my $id (@ids_total) {
			print $OUT_total "$samples_total{$id}{Gender}\t$samples_total{$id}{Size}\n";
			}
		$total_mean_diff = rscript_calc_means($total_outfile,"Female","Male");
		print $OUTI "The CNV data for the individual totals from all available CNVs.\n";
		($info_data_all) = rscript_extract_info($total_outfile,$outfile_final_info);
		#print $LOG1 "The CNV data for the individual totals from all available CNVs.\n".$info_data_all."\n";
		print "\nThere were $count_group1 females CNVs and $count_group2 males CNVs that met criteria for all CNVs.\n";
		print $LOG1 "\nThere were $count_group1 females CNVs and $count_group2 males CNVs that met criteria for all CNVs.\n The mean difference for the groups were $total_mean_diff.\n\nLooking at exonarray available gene CNVs.\n";
		print $OUTF "All CNVs\tComparing Genders of $type\t$total_mean_diff\t$count_group1\t$count_group2\t$totalind1\t$totalind2\n";
		print "Looking at exonarray available gene CNVs.\n"; 
		print "There are ".scalar @exon_keys." CNVs from the exonarray available.\n";
		foreach my $key (@exon_keys){
			$samples_exon{$key}{CNV_ID} = $samples_exon{$key}{$output_headers[$id_col_array]}."_".$samples_exon{$key}{$output_headers[$size_col_array]};
#			print "$key is $samples_exon{$key}{$output_headers[$id_col_array]}.\n";
			if (grep {$samples_exon{$key}{$output_headers[$id_col_array]} eq $_} @excluded) {
				#print "$samples_exon{$key}{$output_headers[$id_col_array]} is excluded.\n";
				next;
				}
			if ($samples_exon{$key}{$output_headers[$type_col1]} eq $type1) {
				if ($samples_exon{$key}{$output_headers[$gender_col_array]} eq $male_gender) {
					push @exon2, $key;
					$count_exon2+=1;
					#print $OUT_temp1 "Male\t".$samples_exon{$key}{$output_headers[$size_col_array]}."\n";
					#print "Male\t$samples_exon{$key}{$output_headers[$id_col_array]}\n";
					if (grep {$samples_exon{$key}{$output_headers[$id_col_array]} eq $_} @ids_exon_total) {
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Size} = $samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Size}+$samples_exon{$key}{$output_headers[$size_col_array]};
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Genes} = $samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Genes}.",".$samples_exon{$key}{$output_headers[$gene_col_array]};
						}
					else {
						push @ids_exon_total, $samples_exon{$key}{$output_headers[$id_col_array]};
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{ID} = $samples_exon{$key}{$output_headers[$id_col_array]};
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Pheno} = $samples_exon{$key}{$output_headers[$aff_col_array]};
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Type} = $samples_exon{$key}{$output_headers[$type_col_array]};
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Gender} = $samples_exon{$key}{$output_headers[$gender_col_array]};
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Size} = $samples_exon{$key}{$output_headers[$size_col_array]};
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Genes} = $samples_exon{$key}{$output_headers[$gene_col_array]};
						$totalexonind2+=1;
						}
					}
				if ($samples_exon{$key}{$output_headers[$gender_col_array]} eq $female_gender) {
					push @exon1, $key; 
					$count_exon1+=1;
					#print $OUT_temp1 "Female\t".$samples_exon{$key}{$output_headers[$size_col_array]}."\n";
					#print "Female\t$samples_exon{$key}{$output_headers[$id_col_array]}\n";
					if (grep {$samples_exon{$key}{$output_headers[$id_col_array]} eq $_} @ids_exon_total) {
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Size} = $samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Size}+$samples_exon{$key}{$output_headers[$size_col_array]};
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Genes} = $samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Genes}.",".$samples_exon{$key}{$output_headers[$gene_col_array]};
						}
					else {
						push @ids_exon_total, $samples_exon{$key}{$output_headers[$id_col_array]};
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{ID} = $samples_exon{$key}{$output_headers[$id_col_array]};
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Pheno} = $samples_exon{$key}{$output_headers[$aff_col_array]};
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Type} = $samples_exon{$key}{$output_headers[$type_col_array]};
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Gender} = $samples_exon{$key}{$output_headers[$gender_col_array]};
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Size} = $samples_exon{$key}{$output_headers[$size_col_array]};
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Genes} = $samples_exon{$key}{$output_headers[$gene_col_array]};
						$totalexonind1+=1;
						}
					}
				}
			}
		foreach my $id (@ids_exon_total) {
			print $OUT_total1 "$samples_exon_total{$id}{Gender}\t$samples_exon_total{$id}{Size}\n";
			}
		$total_exon_mean_diff = rscript_calc_means($total_exon_outfile,"Female","Male");
		print $OUTI "The CNV data for the individual totals from all exonarray available CNVs.\n";
		($info_data_exon) = rscript_extract_info($total_exon_outfile,$outfile_final_info);
		#print $LOG1 "The CNV data for the individual totals from all exonarray available CNVs.\n".$info_data_exon."\n";
		print "There were $count_exon1 females CNVs and $count_exon2 males CNVs that met criteria for exonarray CNVs.\n";
		print $LOG1 "There were $count_exon1 females CNVs and $count_exon2 males CNVs that met criteria for exonarray CNVs.\nThe mean difference for the exonarray available CNVs is $total_exon_mean_diff.\n\nLooking at criteria positive gene CNVs.\n";
		print $OUTF "Exonarray CNVs\tComparing Genders of $type\t$total_exon_mean_diff\t$count_exon1\t$count_exon2\t$totalexonind1\t$totalexonind2\n";
		print "Looking at criteria positive gene CNVs.\n";
		foreach my $key (@criteria_keys){
			$samples_criteria{$key}{CNV_ID} = $samples_criteria{$key}{$output_headers[$id_col_array]}."_".$samples_criteria{$key}{$output_headers[$size_col_array]};
			if (grep {$samples_criteria{$key}{$output_headers[$id_col_array]} eq $_} @excluded) {
				next;
				} 
			if ($samples_criteria{$key}{$output_headers[$type_col1]} eq $type1) {
				if ($samples_criteria{$key}{$output_headers[$gender_col_array]} eq $male_gender) {
				push @criteria2, $key;
				$count_crit2+=1;
				push @criteria_group, $key;
				$criteria_count+=1;
				push @check_crit_cnv_id, $samples_criteria{$key}{CNV_ID};
				#print $OUT_temp2 "Male\t".$samples_criteria{$key}{$output_headers[$size_col_array]}."\n";
					if (grep {$samples_criteria{$key}{$output_headers[$id_col_array]} eq $_} @ids_crit_total) {
						#print "$samples{$key}{$output_headers[$id_col_array]}\t";
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Size} = $samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Size}+$samples_criteria{$key}{$output_headers[$size_col_array]};
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Genes} = $samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Genes}.",".$samples_criteria{$key}{$output_headers[$gene_col_array]};
						$numcnvs2+=1;
						#print "$numcnvs2\n";
						}
					else {
						#print "2\t$samples{$key}{$output_headers[$id_col_array]}\t";
						push @ids_crit_total, $samples_criteria{$key}{$output_headers[$id_col_array]};
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{ID} = $samples_criteria{$key}{$output_headers[$id_col_array]};
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Pheno} = $samples_criteria{$key}{$output_headers[$aff_col_array]};
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Type} = $samples_criteria{$key}{$output_headers[$type_col_array]};
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Gender} = $samples_criteria{$key}{$output_headers[$gender_col_array]};
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Size} = $samples_criteria{$key}{$output_headers[$size_col_array]};
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Genes} = $samples_criteria{$key}{$output_headers[$gene_col_array]};
						$numcnvs2+=1;
						$totalcritind2+=1;
						#print "$numcnvs2\t$totalcritind2\n"
						}
					}
				if ($samples_criteria{$key}{$output_headers[$gender_col_array]} eq $female_gender) {
					push @criteria1, $key; 
					$count_crit1+=1;
					push @criteria_group, $key;
					$criteria_count+=1;
					push @check_crit_cnv_id, $samples_criteria{$key}{CNV_ID};
				#	print $OUT_temp2 "Female\t".$samples_criteria{$key}{$output_headers[$size_col_array]}."\n";
					if (grep {$samples_criteria{$key}{$output_headers[$id_col_array]} eq $_} @ids_crit_total) {
						#print "$samples{$key}{$output_headers[$id_col_array]}\t";
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Size} = $samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Size}+$samples_criteria{$key}{$output_headers[$size_col_array]};
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Genes} = $samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Genes}.",".$samples_criteria{$key}{$output_headers[$gene_col_array]};
						$numcnvs1+=1;
						#print "$numcnvs1\n"
						}
					else {
						#print "1\t$samples{$key}{$output_headers[$id_col_array]}\t";
						push @ids_crit_total, $samples_criteria{$key}{$output_headers[$id_col_array]};
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{ID} = $samples_criteria{$key}{$output_headers[$id_col_array]};
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Pheno} = $samples_criteria{$key}{$output_headers[$aff_col_array]};
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Type} = $samples_criteria{$key}{$output_headers[$type_col_array]};
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Gender} = $samples_criteria{$key}{$output_headers[$gender_col_array]};
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Size} = $samples_criteria{$key}{$output_headers[$size_col_array]};
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Genes} = $samples_criteria{$key}{$output_headers[$gene_col_array]};
						$numcnvs1+=1;
						$totalcritind1+=1;
						#print "$numcnvs1\t$totalcritind1\n"
						}
					}
				}
			}
		foreach my $id (@ids_crit_total) {
			print $OUT_total2 "$samples_crit_total{$id}{Gender}\t$samples_crit_total{$id}{Size}\n";
			}
		$total_crit_mean_diff = rscript_calc_means($total_crit_outfile,"Female","Male");
		print $OUTI "The CNV data for the individual totals from all criteria available CNVs.\n";
		($info_data_criteria) = rscript_extract_info($total_crit_outfile,$outfile_final_info);
		print "The numcnvs1 is $numcnvs1 and numcnvs2 is $numcnvs2.\n";
		print "There were $count_crit1 females CNVs and $count_crit2 males CNVs that met criteria for Criteria CNVs.\n";
		print "There are $criteria_count overall samples to use for random selection in permutations.\n";
		print $LOG1 "There were $count_crit1 females CNVs and $count_crit2 males CNVs that met criteria for Criteria CNVs.\nThere are $criteria_count overall samples to use for random selection in permutations.\nThe mean difference in comparing the criteria met CNVs per individual is $total_crit_mean_diff.\n\n";
		print $OUTF "Criteria CNVs\tComparing Genders of $type\t$total_crit_mean_diff\t$count_crit1\t$count_crit2\t$totalcritind1\t$totalcritind2\n";
		print $OUT "0\tCriteria CNVs\tComparing Genders of $type\t$total_crit_mean_diff\n";
		print $OUTC "0\tCriteria CNVs\tComparing Genders of $type\t$total_crit_mean_diff\n";
		print $OUTE "0\tCriteria CNVs\tComparing Genders of $type\t$total_crit_mean_diff\n";
		print $OUTEG "0\tCriteria CNVs\tComparing Genders of $type\t$total_crit_mean_diff\n";
		print $OUTG "0\tCriteria CNVs\tComparing Genders of $type\t$total_crit_mean_diff\n";
		$comparison_text = "Comparing Genders of $type";
		if ($number_for_groups) {	
			print $OUTN "0\tCriteria CNVs $num_gp1 vs $num_gp2\tComparing Genders of $type\t$numbers_mean\n";
			print $OUTF "Criteria CNVs $num_gp1 vs $num_gp2\tComparing Genders of $type\t$numbers_mean\n";
			print "Looking at $num_gp1 females CNVs and $num_gp2 males CNVs for comparison to mean $numbers_mean.\n";
			}
	
		print $OUTFT "Sample_ID\tPhenotype\tSample_Type\tGender\tCNV_Size\tGenes\n";
		foreach my $ids (@ids_crit_total){
			print $OUTFT "$samples_crit_total{$ids}{ID}\t$samples_crit_total{$ids}{Pheno}\t$samples_crit_total{$ids}{Type}\t$samples_crit_total{$ids}{Gender}\t$samples_crit_total{$ids}{Size}\t$samples_crit_total{$ids}{Genes}\n";
			}
		print $LOG1 "The CNV_IDs for the criteria region is @check_crit_cnv_id\n";
		}
		
	elsif ($sex && $affected && !$gender) {
		print "Comparing the affected status for $sex for $type1 and $type2.\n";
		print $OUT_temp "Type\tSize\n";
		print $OUT_temp1 "Type\tSize\n";
		print $OUT_temp2 "Type\tSize\n";
		print $OUT_total "Type\tSize\n";
		print $OUT_total1 "Type\tSize\n";
		print $OUT_total2 "Type\tSize\n";
		foreach my $key (@sample_keys){
			#print "$samples{$key}{$output_headers[$id_col_array]}.\n";
			if (grep {$samples{$key}{$output_headers[$id_col_array]} eq $_} @excluded) {next;}  
			if ($samples{$key}{$output_headers[$gender_col_array]} eq $sex) {
				#print "$samples{$key}{$output_headers[$id_col_array]} is $sex.\t$samples{$key}{$output_headers[$aff_col_array]}\n";
				if ($samples{$key}{$output_headers[$aff_col_array]} eq $type1) {
					push @group1, $key;
					$count_group1+=1;
					#print "$samples{$key}{$output_headers[$id_col_array]} is $type1.\n";
					print $OUT_temp $type1."\t".$samples{$key}{$output_headers[$size_col_array]}."\n";
					if (grep {$samples{$key}{$output_headers[$id_col_array]} eq $_} @ids_total) {
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Size} = $samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Size}+$samples{$key}{$output_headers[$size_col_array]};
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Genes} = $samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Genes}.",".$samples{$key}{$output_headers[$gene_col_array]};
						}
					else {
						#print "First time seen $samples{$key}{$output_headers[$id_col_array]}\n";
						push @ids_total, $samples{$key}{$output_headers[$id_col_array]};
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{ID} = $samples{$key}{$output_headers[$id_col_array]};
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Pheno} = $samples{$key}{$output_headers[$aff_col_array]};
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Type} = $samples{$key}{$output_headers[$type_col_array]};
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Gender} = $samples{$key}{$output_headers[$gender_col_array]};
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Size} = $samples{$key}{$output_headers[$size_col_array]};
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Genes} = $samples{$key}{$output_headers[$gene_col_array]};
						$totalind1+=1;
						}
				}
				if ($samples{$key}{$output_headers[$aff_col_array]} eq $type2) {
					push @group2, $key; 
					$count_group2+=1;
					#print "$samples{$key}{$output_headers[$id_col_array]} is $type1.\n";
					print $OUT_temp $type2."\t".$samples{$key}{$output_headers[$size_col_array]}."\n";
					if (grep {$samples{$key}{$output_headers[$id_col_array]} eq $_} @ids_total) {
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Size} = $samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Size}+$samples{$key}{$output_headers[$size_col_array]};
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Genes} = $samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Genes}.",".$samples{$key}{$output_headers[$gene_col_array]};
						}
					else {
						#print "First time seen $samples{$key}{$output_headers[$id_col_array]}\n";
						push @ids_total, $samples{$key}{$output_headers[$id_col_array]};
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{ID} = $samples{$key}{$output_headers[$id_col_array]};
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Pheno} = $samples{$key}{$output_headers[$aff_col_array]};
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Type} = $samples{$key}{$output_headers[$type_col_array]};
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Gender} = $samples{$key}{$output_headers[$gender_col_array]};
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Size} = $samples{$key}{$output_headers[$size_col_array]};
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Genes} = $samples{$key}{$output_headers[$gene_col_array]};
						$totalind2+=1;
						}
					}
				}
			}
		foreach my $id (@ids_total) {
			print $OUT_total "$samples_total{$id}{Pheno}\t$samples_total{$id}{Size}\n";
			}
		$total_mean_diff = rscript_calc_means($total_outfile,$type1,$type2);
		print $OUTI "The CNV data for the individual totals from all available CNVs.\n";
		($info_data_all) = rscript_extract_info($total_outfile,$outfile_final_info);
		#print $LOG1 "The CNV data for the individual totals from all available CNVs.\n$info_data_all\n";
		print $OUTF "All CNVs\tComparing Affected Status of $sex\t$total_mean_diff\t$count_group1\t$count_group2\t$totalind1\t$totalind2\n";
		foreach my $key (@exon_keys){
			$samples_exon{$key}{CNV_ID} = $samples_exon{$key}{$output_headers[$id_col_array]}."_".$samples_exon{$key}{$output_headers[$size_col_array]};
			if (grep {$samples_exon{$key}{$output_headers[$id_col_array]} eq $_} @excluded) {next;}
			if ($samples_exon{$key}{$output_headers[$gender_col_array]} eq $sex) {
				if ($samples_exon{$key}{$output_headers[$aff_col_array]} eq $type1) {
					push @exon1, $key;
					$count_exon1+=1;
					print $OUT_temp1 $type1."\t".$samples_exon{$key}{$output_headers[$size_col_array]}."\n";
					if (grep {$samples_exon{$key}{$output_headers[$id_col_array]} eq $_} @ids_exon_total) {
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Size} = $samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Size}+$samples_exon{$key}{$output_headers[$size_col_array]};
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Genes} = $samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Genes}.",".$samples_exon{$key}{$output_headers[$gene_col_array]};
						}
					else {
						push @ids_exon_total, $samples_exon{$key}{$output_headers[$id_col_array]};
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{ID} = $samples_exon{$key}{$output_headers[$id_col_array]};
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Pheno} = $samples_exon{$key}{$output_headers[$aff_col_array]};
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Type} = $samples_exon{$key}{$output_headers[$type_col_array]};
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Gender} = $samples_exon{$key}{$output_headers[$gender_col_array]};
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Size} = $samples_exon{$key}{$output_headers[$size_col_array]};
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Genes} = $samples_exon{$key}{$output_headers[$gene_col_array]};
						$totalexonind1+=1;
						}
					}
				if ($samples_exon{$key}{$output_headers[$aff_col_array]} eq $type2) {
					push @exon2, $key; 
					$count_exon2+=1;
					print $OUT_temp1 $type2."\t".$samples_exon{$key}{$output_headers[$size_col_array]}."\n";
					if (grep {$samples_exon{$key}{$output_headers[$id_col_array]} eq $_} @ids_exon_total) {
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Size} = $samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Size}+$samples_exon{$key}{$output_headers[$size_col_array]};
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Genes} = $samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Genes}.",".$samples_exon{$key}{$output_headers[$gene_col_array]};
						}
					else {
						push @ids_exon_total, $samples_exon{$key}{$output_headers[$id_col_array]};
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{ID} = $samples_exon{$key}{$output_headers[$id_col_array]};
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Pheno} = $samples_exon{$key}{$output_headers[$aff_col_array]};
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Type} = $samples_exon{$key}{$output_headers[$type_col_array]};
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Gender} = $samples_exon{$key}{$output_headers[$gender_col_array]};
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Size} = $samples_exon{$key}{$output_headers[$size_col_array]};
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Genes} = $samples_exon{$key}{$output_headers[$gene_col_array]};
						$totalexonind2+=1;
						}
					}
				}
			}
		foreach my $id (@ids_exon_total) {
			print $OUT_total1 "$samples_exon_total{$id}{Pheno}\t$samples_exon_total{$id}{Size}\n";
			}
		$total_exon_mean_diff = rscript_calc_means($total_exon_outfile,$type1,$type2);
		print $OUTI "The CNV data for the individual totals from all exonarray available CNVs.\n";
		($info_data_exon) = rscript_extract_info($total_exon_outfile,$outfile_final_info);
		#print $LOG1 "The CNV data for the individual totals from all exonarray available CNVs.\n$info_data_exon\n";
		print $OUTF "Exonarray CNVs\tComparing Affected Status of $sex\t$total_exon_mean_diff\t$count_exon1\t$count_exon2\t$totalexonind1\t$totalexonind2\n";
		foreach my $key (@criteria_keys){
			$samples_criteria{$key}{CNV_ID} = $samples_criteria{$key}{$output_headers[$id_col_array]}."_".$samples_criteria{$key}{$output_headers[$size_col_array]};
			if (grep {$samples_criteria{$key}{$output_headers[$id_col_array]} eq $_} @excluded) {next;} 
			if ($samples_criteria{$key}{$output_headers[$gender_col_array]} eq $sex) {
				if ($samples_criteria{$key}{$output_headers[$aff_col_array]} eq $type1) {
					push @criteria1, $key;
					$count_crit1+=1;
					push @criteria_group, $key;
					$criteria_count+=1;
					push @check_crit_cnv_id, $samples_criteria{$key}{CNV_ID};
					print $OUT_temp2 $type1."\t".$samples_criteria{$key}{$output_headers[$size_col_array]}."\n";
					if (grep {$samples_criteria{$key}{$output_headers[$id_col_array]} eq $_} @ids_crit_total) {
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Size} = $samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Size}+$samples_criteria{$key}{$output_headers[$size_col_array]};
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Genes} = $samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Genes}.",".$samples_criteria{$key}{$output_headers[$gene_col_array]};
						$numcnvs1+=1;
						}
					else {
						push @ids_crit_total, $samples_criteria{$key}{$output_headers[$id_col_array]};
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{ID} = $samples_criteria{$key}{$output_headers[$id_col_array]};
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Pheno} = $samples_criteria{$key}{$output_headers[$aff_col_array]};
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Type} = $samples_criteria{$key}{$output_headers[$type_col_array]};
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Gender} = $samples_criteria{$key}{$output_headers[$gender_col_array]};
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Size} = $samples_criteria{$key}{$output_headers[$size_col_array]};
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Genes} = $samples_criteria{$key}{$output_headers[$gene_col_array]};
						$numcnvs1+=1;
						$totalcritind1+=1;
						}
					}
				if ($samples_criteria{$key}{$output_headers[$aff_col_array]} eq $type2) {
					push @criteria2, $key; 
					$count_crit2+=1;
					push @criteria_group, $key;
					$criteria_count+=1;
					push @check_crit_cnv_id, $samples_criteria{$key}{CNV_ID};
					print $OUT_temp2 $type2."\t".$samples_criteria{$key}{$output_headers[$size_col_array]}."\n";
					if (grep {$samples_criteria{$key}{$output_headers[$id_col_array]} eq $_} @ids_crit_total) {
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Size} = $samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Size}+$samples_criteria{$key}{$output_headers[$size_col_array]};
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Genes} = $samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Genes}.",".$samples_criteria{$key}{$output_headers[$gene_col_array]};
						$numcnvs2+=1;
						}
					else {
						push @ids_crit_total, $samples_criteria{$key}{$output_headers[$id_col_array]};
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{ID} = $samples_criteria{$key}{$output_headers[$id_col_array]};
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Pheno} = $samples_criteria{$key}{$output_headers[$aff_col_array]};
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Type} = $samples_criteria{$key}{$output_headers[$type_col_array]};
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Gender} = $samples_criteria{$key}{$output_headers[$gender_col_array]};
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Size} = $samples_criteria{$key}{$output_headers[$size_col_array]};
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Genes} = $samples_criteria{$key}{$output_headers[$gene_col_array]};
						$numcnvs2+=1;
						$totalcritind2+=1;
						}
					}
				}
			}
		foreach my $id (@ids_crit_total) {
			print $OUT_total2 "$samples_crit_total{$id}{Pheno}\t$samples_crit_total{$id}{Size}\n";
			}		
		$total_crit_mean_diff = rscript_calc_means($total_crit_outfile,$type1,$type2);
		print $OUTI "The CNV data for the individual totals from criteria available CNVs.\n";
		($info_data_criteria) = rscript_extract_info($total_crit_outfile,$outfile_final_info);
		print $OUTF "Criteria CNVs\tComparing Affected Status of $sex\t$total_crit_mean_diff\t$count_crit1\t$count_crit2\t$totalcritind1\t$totalcritind2\n";
		print $OUT "0\tCriteria CNVs\tComparing Affected Status of $sex\t$total_crit_mean_diff\n";
		print $OUTC "0\tCriteria CNVs\tComparing Affected Status of $sex\t$total_crit_mean_diff\n";
		print $OUTE "0\tCriteria CNVs\tComparing Affected Status of $sex\t$total_crit_mean_diff\n";
		print $OUTEG "0\tCriteria CNVs\tComparing Affected Status of $sex\t$total_crit_mean_diff\n";
		print $OUTG "0\tCriteria CNVs\tComparing Affected Status of $sex\t$total_crit_mean_diff\n";
		$comparison_text = "Comparing Affected Status of $sex";
		if ($number_for_groups) {	
			print $OUTN "0\tCriteria CNVs $num_gp1 vs $num_gp2\tComparing Affected Status of $sex\t$numbers_mean\n";
			print $OUTF "Criteria CNVs $num_gp1 vs $num_gp2\tComparing Affected Status of $sex\t$numbers_mean\n";}
		
		print $OUTFT "Sample_ID\tPhenotype\tSample_Type\tGender\tCNV_Size\tGenes\n";
		foreach my $ids (@ids_crit_total){
			print $OUTFT "$samples_crit_total{$ids}{ID}\t$samples_crit_total{$ids}{Pheno}\t$samples_crit_total{$ids}{Type}\t$samples_crit_total{$ids}{Gender}\t$samples_crit_total{$ids}{Size}\t$samples_crit_total{$ids}{Genes}\n";
			}
		}
		
	elsif ($sex && $type && !$gender) {
		print "Comparing $type for $sex.\n";
		print $OUT_temp "Type\tSize\n";
		print $OUT_temp1 "Type\tSize\n";
		print $OUT_temp2 "Type\tSize\n";
		print $OUT_total "Type\tSize\n";
		print $OUT_total1 "Type\tSize\n";
		print $OUT_total2 "Type\tSize\n";
		foreach my $key (@sample_keys){
			if (grep {$samples{$key}{$output_headers[$id_col_array]} eq $_} @excluded) {next;} 
			if ($samples{$key}{$output_headers[$gender_col_array]} =~ /$sex/i) {
				if ($samples{$key}{$output_headers[$type_col1]} =~ /$type1/i) {
					push @group1, $key;
					$count_group1+=1;
					print $OUT_temp $type1."\t".$samples{$key}{$output_headers[$size_col_array]}."\n";
					if (grep {$samples{$key}{$output_headers[$id_col_array]} eq $_} @ids_total) {
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Size} = $samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Size}+$samples{$key}{$output_headers[$size_col_array]};
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Genes} = $samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Genes}.",".$samples{$key}{$output_headers[$gene_col_array]};
						}
					else {
						push @ids_total, $samples{$key}{$output_headers[$id_col_array]};
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{ID} = $samples{$key}{$output_headers[$id_col_array]};
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Pheno} = $samples{$key}{$output_headers[$aff_col_array]};
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Type} = $samples{$key}{$output_headers[$type_col_array]};
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Gender} = $samples{$key}{$output_headers[$gender_col_array]};
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Size} = $samples{$key}{$output_headers[$size_col_array]};
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Genes} = $samples{$key}{$output_headers[$gene_col_array]};
						$totalind1+=1;
						}
					}
				if ($samples{$key}{$output_headers[$type_col2]} =~ /$type2/i) {
					push @group2, $key; 
					$count_group2+=1;
					print $OUT_temp $type2."\t".$samples{$key}{$output_headers[$size_col_array]}."\n";
					if (grep {$samples{$key}{$output_headers[$id_col_array]} eq $_} @ids_total) {
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Size} = $samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Size}+$samples{$key}{$output_headers[$size_col_array]};
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Genes} = $samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Genes}.",".$samples{$key}{$output_headers[$gene_col_array]};
						}
					else {
						push @ids_total, $samples{$key}{$output_headers[$id_col_array]};
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{ID} = $samples{$key}{$output_headers[$id_col_array]};
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Pheno} = $samples{$key}{$output_headers[$aff_col_array]};
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Type} = $samples{$key}{$output_headers[$type_col_array]};
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Gender} = $samples{$key}{$output_headers[$gender_col_array]};
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Size} = $samples{$key}{$output_headers[$size_col_array]};
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Genes} = $samples{$key}{$output_headers[$gene_col_array]};
						$totalind2+=1;
						}
					}
				}
			}
		foreach my $id (@ids_total) {
			print $OUT_total "$samples_total{$id}{Type}\t$samples_total{$id}{Size}\n";
			}
		$total_mean_diff = rscript_calc_means($total_outfile,$type1,$type2);
		print $OUTI "The CNV data for the individual totals from all available CNVs.\n";
		($info_data_all) = rscript_extract_info($total_outfile,$outfile_final_info);
		#print $LOG1 "The CNV data for the individual totals from all available CNVs.\n$info_data_all\n";
		print $OUTF "All CNVs\tComparing $type1 and $type2 of $sex\t$total_mean_diff\t$count_group1\t$count_group2\t$totalind1\t$totalind2\n";
		foreach my $key (@exon_keys){
			$samples_exon{$key}{CNV_ID} = $samples_exon{$key}{$output_headers[$id_col_array]}."_".$samples_exon{$key}{$output_headers[$size_col_array]};
			if (grep {$samples_exon{$key}{$output_headers[$id_col_array]} eq $_} @excluded) {next;} 
			if ($samples_exon{$key}{$output_headers[$gender_col_array]} =~ /$sex/i) {
				if ($samples_exon{$key}{$output_headers[$type_col1]} =~ /$type1/i) {
					push @exon1, $key;
					$count_exon1+=1;
					print $OUT_temp1 $type1."\t".$samples_exon{$key}{$output_headers[$size_col_array]}."\n";
					if (grep {$samples_exon{$key}{$output_headers[$id_col_array]} eq $_} @ids_exon_total) {
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Size} = $samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Size}+$samples_exon{$key}{$output_headers[$size_col_array]};
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Genes} = $samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Genes}.",".$samples_exon{$key}{$output_headers[$gene_col_array]};
						}
					else {
						push @ids_exon_total, $samples_exon{$key}{$output_headers[$id_col_array]};
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{ID} = $samples_exon{$key}{$output_headers[$id_col_array]};
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Pheno} = $samples_exon{$key}{$output_headers[$aff_col_array]};
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Type} = $samples_exon{$key}{$output_headers[$type_col_array]};
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Gender} = $samples_exon{$key}{$output_headers[$gender_col_array]};
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Size} = $samples_exon{$key}{$output_headers[$size_col_array]};
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Genes} = $samples_exon{$key}{$output_headers[$gene_col_array]};
						$totalexonind1+=1;
						}
					}
				if ($samples_exon{$key}{$output_headers[$type_col2]} =~ /$type2/i) {
					push @exon2, $key; 
					$count_exon2+=1;
					print $OUT_temp1 $type2."\t".$samples_exon{$key}{$output_headers[$size_col_array]}."\n";
					if (grep {$samples_exon{$key}{$output_headers[$id_col_array]} eq $_} @ids_exon_total) {
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Size} = $samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Size}+$samples_exon{$key}{$output_headers[$size_col_array]};
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Genes} = $samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Genes}.",".$samples_exon{$key}{$output_headers[$gene_col_array]};
						}
					else {
						push @ids_exon_total, $samples_exon{$key}{$output_headers[$id_col_array]};
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{ID} = $samples_exon{$key}{$output_headers[$id_col_array]};
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Pheno} = $samples_exon{$key}{$output_headers[$aff_col_array]};
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Type} = $samples_exon{$key}{$output_headers[$type_col_array]};
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Gender} = $samples_exon{$key}{$output_headers[$gender_col_array]};
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Size} = $samples_exon{$key}{$output_headers[$size_col_array]};
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Genes} = $samples_exon{$key}{$output_headers[$gene_col_array]};
						$totalexonind2+=1;
						}
					}
				}
			}
		foreach my $id (@ids_exon_total) {
			print $OUT_total1 "$samples_exon_total{$id}{Type}\t$samples_exon_total{$id}{Size}\n";
			}
		$total_exon_mean_diff = rscript_calc_means($total_exon_outfile,$type1,$type2);
		print $OUTI "The CNV data for the individual totals from all exonarray available CNVs.\n";
		($info_data_exon) = rscript_extract_info($total_exon_outfile,$outfile_final_info);
		print $OUTF "Exonarray CNVs\tComparing $type1 and $type2 of $sex\t$total_exon_mean_diff\t$count_exon1\t$count_exon2\t$totalexonind1\t$totalexonind2\n";
		foreach my $key (@criteria_keys){
			$samples_criteria{$key}{CNV_ID} = $samples_criteria{$key}{$output_headers[$id_col_array]}."_".$samples_criteria{$key}{$output_headers[$size_col_array]};
			if (grep {$samples_criteria{$key}{$output_headers[$id_col_array]} eq $_} @excluded) {next;} 
			if ($samples_criteria{$key}{$output_headers[$gender_col_array]} =~ /$sex/i) {
				if ($samples_criteria{$key}{$output_headers[$type_col1]} =~ /$type1/i) {
					push @criteria1, $key;
					$count_crit1+=1;
					push @criteria_group, $key;
					$criteria_count+=1;
					push @check_crit_cnv_id, $samples_criteria{$key}{CNV_ID};
					print $OUT_temp2 $type1."\t".$samples_criteria{$key}{$output_headers[$size_col_array]}."\n";
					if (grep {$samples_criteria{$key}{$output_headers[$id_col_array]} eq $_} @ids_crit_total) {
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Size} = $samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Size}+$samples_criteria{$key}{$output_headers[$size_col_array]};
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Genes} = $samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Genes}.",".$samples_criteria{$key}{$output_headers[$gene_col_array]};
						$numcnvs1+=1;
						}
					else {
						push @ids_crit_total, $samples_criteria{$key}{$output_headers[$id_col_array]};
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{ID} = $samples_criteria{$key}{$output_headers[$id_col_array]};
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Pheno} = $samples_criteria{$key}{$output_headers[$aff_col_array]};
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Type} = $samples_criteria{$key}{$output_headers[$type_col_array]};
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Gender} = $samples_criteria{$key}{$output_headers[$gender_col_array]};
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Size} = $samples_criteria{$key}{$output_headers[$size_col_array]};
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Genes} = $samples_criteria{$key}{$output_headers[$gene_col_array]};
						$numcnvs1+=1;
						$totalcritind1+=1;
						}
					}
				if ($samples_criteria{$key}{$output_headers[$type_col2]} =~ /$type2/i) {
					push @criteria2, $key; 
					$count_crit2+=1;
					push @criteria_group, $key;
					$criteria_count+=1;
					push @check_crit_cnv_id, $samples_criteria{$key}{CNV_ID};
					print $OUT_temp2 $type2."\t".$samples_criteria{$key}{$output_headers[$size_col_array]}."\n";
					if (grep {$samples_criteria{$key}{$output_headers[$id_col_array]} eq $_} @ids_crit_total) {
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Size} = $samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Size}+$samples_criteria{$key}{$output_headers[$size_col_array]};
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Genes} = $samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Genes}.",".$samples_criteria{$key}{$output_headers[$gene_col_array]};
						$numcnvs2+=1;
						}
					else {
						push @ids_crit_total, $samples_criteria{$key}{$output_headers[$id_col_array]};
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{ID} = $samples_criteria{$key}{$output_headers[$id_col_array]};
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Pheno} = $samples_criteria{$key}{$output_headers[$aff_col_array]};
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Type} = $samples_criteria{$key}{$output_headers[$type_col_array]};
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Gender} = $samples_criteria{$key}{$output_headers[$gender_col_array]};
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Size} = $samples_criteria{$key}{$output_headers[$size_col_array]};
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Genes} = $samples_criteria{$key}{$output_headers[$gene_col_array]};
						$numcnvs2+=1;
						$totalcritind2+=1;
						}
					}
				}
			}
		foreach my $id (@ids_crit_total) {
			print $OUT_total2 "$samples_crit_total{$id}{Type}\t$samples_crit_total{$id}{Size}\n";
			}	
		$total_crit_mean_diff = rscript_calc_means($total_crit_outfile,$type1,$type2);
		print $OUTI "The CNV data for the individual totals from all criteria available CNVs.\n";
		($info_data_criteria) = rscript_extract_info($total_crit_outfile,$outfile_final_info);
		print $OUTF "Criteria CNVs\tComparing $type1 and $type2 of $sex\t$total_crit_mean_diff\t$count_crit1\t$count_crit2\t$totalcritind1\t$totalcritind2\n";
		print $OUT "0\tCriteria CNVs\tComparing $type1 and $type2 of $sex\t$total_crit_mean_diff\n";
		print $OUTC "0\tCriteria CNVs\tComparing $type1 and $type2 of $sex\t$total_crit_mean_diff\n";
		print $OUTE "0\tCriteria CNVs\tComparing $type1 and $type2 of $sex\t$total_crit_mean_diff\n";
		print $OUTEG "0\tCriteria CNVs\tComparing $type1 and $type2 of $sex\t$total_crit_mean_diff\n";
		print $OUTG "0\tCriteria CNVs\tComparing $type1 and $type2 of $sex\t$total_crit_mean_diff\n";
		$comparison_text = "Comparing $type1 and $type2";
		if ($number_for_groups) {	
			print $OUTN "0\tCriteria CNVs $num_gp1 vs $num_gp2\tComparing $type1 and $type2 of $sex\t$numbers_mean\n";
			print $OUTF "Criteria CNVs $num_gp1 vs $num_gp2\tComparing $type1 and $type2 of $sex\t$numbers_mean\n";}
		print $OUTFT "Sample_ID\tPhenotype\tSample_Type\tGender\tCNV_Size\tGenes\n";
		foreach my $ids (@ids_crit_total){
			print $OUTFT "$samples_crit_total{$ids}{ID}\t$samples_crit_total{$ids}{Pheno}\t$samples_crit_total{$ids}{Type}\t$samples_crit_total{$ids}{Gender}\t$samples_crit_total{$ids}{Size}\t$samples_crit_total{$ids}{Genes}\n";
			}
		}
		
	elsif ($affected) {
		print "Comparing the affected status of the samples.\n";
		print $OUT_temp "Type\tSize\n";
		print $OUT_temp1 "Type\tSize\n";
		print $OUT_temp2 "Type\tSize\n";
		print $OUT_total "Type\tSize\n";
		print $OUT_total1 "Type\tSize\n";
		print $OUT_total2 "Type\tSize\n";
		foreach my $key (@sample_keys){
			if (grep {$samples{$key}{$output_headers[$id_col_array]} eq $_} @excluded) {print "This sample is excluded.\n";next;}  
			if ($samples{$key}{$output_headers[$aff_col_array]} =~ /$type1/i) {
				push @group1, $key;
				$count_group1+=1;
				print $OUT_temp $type1."\t".$samples{$key}{$output_headers[$size_col_array]}."\n";
				if (grep {$samples{$key}{$output_headers[$id_col_array]} eq $_} @ids_total) {
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Size} = $samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Size}+$samples{$key}{$output_headers[$size_col_array]};
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Genes} = $samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Genes}.",".$samples{$key}{$output_headers[$gene_col_array]};
						}
					else {
						push @ids_total, $samples{$key}{$output_headers[$id_col_array]};
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{ID} = $samples{$key}{$output_headers[$id_col_array]};
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Pheno} = $samples{$key}{$output_headers[$aff_col_array]};
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Type} = $samples{$key}{$output_headers[$type_col_array]};
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Gender} = $samples{$key}{$output_headers[$gender_col_array]};
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Size} = $samples{$key}{$output_headers[$size_col_array]};
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Genes} = $samples{$key}{$output_headers[$gene_col_array]};
						$totalind1+=1;
						}
					}
			if ($samples{$key}{$output_headers[$aff_col_array]} =~ /$type2/i) {
				push @group2, $key; 
				$count_group2+=1;
				print $OUT_temp $type2."\t".$samples{$key}{$output_headers[$size_col_array]}."\n";
				if (grep {$samples{$key}{$output_headers[$id_col_array]} eq $_} @ids_total) {
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Size} = $samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Size}+$samples{$key}{$output_headers[$size_col_array]};
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Genes} = $samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Genes}.",".$samples{$key}{$output_headers[$gene_col_array]};
						}
					else {
						push @ids_total, $samples{$key}{$output_headers[$id_col_array]};
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{ID} = $samples{$key}{$output_headers[$id_col_array]};
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Pheno} = $samples{$key}{$output_headers[$aff_col_array]};
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Type} = $samples{$key}{$output_headers[$type_col_array]};
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Gender} = $samples{$key}{$output_headers[$gender_col_array]};
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Size} = $samples{$key}{$output_headers[$size_col_array]};
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Genes} = $samples{$key}{$output_headers[$gene_col_array]};
						$totalind2+=1;
						}
					}
			}
		foreach my $id (@ids_total) {
			print $OUT_total "$samples_total{$id}{Pheno}\t$samples_total{$id}{Size}\n";
			}
		$total_mean_diff = rscript_calc_means($total_outfile,$type1,$type2);
		print $OUTI "The CNV data for the individual totals for all available CNVs.\n";
		($info_data_all) = rscript_extract_info($total_outfile,$outfile_final_info);
		#print $LOG1 "The CNV data for the individual totals for all available CNVs.\n$info_data_all\n";
		print $OUTF "All CNVs\tComparing $type1 and $type2\t$total_mean_diff\t$count_group1\t$count_group2\t$totalind1\t$totalind2\n";
		foreach my $key (@exon_keys){
			$samples_exon{$key}{CNV_ID} = $samples_exon{$key}{$output_headers[$id_col_array]}."_".$samples_exon{$key}{$output_headers[$size_col_array]};
			if (grep {$samples_exon{$key}{$output_headers[$id_col_array]} eq $_} @excluded) {next;}
			if ($samples_exon{$key}{$output_headers[$aff_col_array]} =~ /$type1/i) {
				push @exon1, $key;
				$count_exon1+=1;
				print $OUT_temp1 $type1."\t".$samples_exon{$key}{$output_headers[$size_col_array]}."\n";
				if (grep {$samples_exon{$key}{$output_headers[$id_col_array]} eq $_} @ids_exon_total) {
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Size} = $samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Size}+$samples_exon{$key}{$output_headers[$size_col_array]};
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Genes} = $samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Genes}.",".$samples_exon{$key}{$output_headers[$gene_col_array]};
						}
					else {
						push @ids_exon_total, $samples_exon{$key}{$output_headers[$id_col_array]};
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{ID} = $samples_exon{$key}{$output_headers[$id_col_array]};
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Pheno} = $samples_exon{$key}{$output_headers[$aff_col_array]};
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Type} = $samples_exon{$key}{$output_headers[$type_col_array]};
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Gender} = $samples_exon{$key}{$output_headers[$gender_col_array]};
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Size} = $samples_exon{$key}{$output_headers[$size_col_array]};
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Genes} = $samples_exon{$key}{$output_headers[$gene_col_array]};
						$totalexonind1+=1;
						}
					}
			if ($samples_exon{$key}{$output_headers[$aff_col_array]} =~ /$type2/i) {
				push @exon2, $key; 
				$count_exon2+=1;
				print $OUT_temp1 $type2."\t".$samples_exon{$key}{$output_headers[$size_col_array]}."\n";
				if (grep {$samples_exon{$key}{$output_headers[$id_col_array]} eq $_} @ids_exon_total) {
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Size} = $samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Size}+$samples_exon{$key}{$output_headers[$size_col_array]};
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Genes} = $samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Genes}.",".$samples_exon{$key}{$output_headers[$gene_col_array]};
						}
					else {
						push @ids_exon_total, $samples_exon{$key}{$output_headers[$id_col_array]};
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{ID} = $samples_exon{$key}{$output_headers[$id_col_array]};
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Pheno} = $samples_exon{$key}{$output_headers[$aff_col_array]};
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Type} = $samples_exon{$key}{$output_headers[$type_col_array]};
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Gender} = $samples_exon{$key}{$output_headers[$gender_col_array]};
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Size} = $samples_exon{$key}{$output_headers[$size_col_array]};
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Genes} = $samples_exon{$key}{$output_headers[$gene_col_array]};
						$totalexonind2+=1;
						}
					}
			}
		foreach my $id (@ids_exon_total) {
			print $OUT_total1 "$samples_exon_total{$id}{Pheno}\t$samples_exon_total{$id}{Size}\n";
			}
		$total_exon_mean_diff = rscript_calc_means($total_exon_outfile,$type1,$type2);
		print $OUTI "The CNV data for the individual totals for all exonarray available CNVs.\n";
		($info_data_exon) = rscript_extract_info($total_exon_outfile,$outfile_final_info);
		#print $LOG1 "The CNV data for the individual totals for all exonarray available CNVs.\n$info_data_exon\n";
		print $OUTF "Exonarray CNVs\tComparing $type1 and $type2\t$total_exon_mean_diff\t$count_exon1\t$count_exon2\t$totalexonind1\t$totalexonind2\n";
		foreach my $key (@criteria_keys){
			$samples_criteria{$key}{CNV_ID} = $samples_criteria{$key}{$output_headers[$id_col_array]}."_".$samples_criteria{$key}{$output_headers[$size_col_array]};
			if (grep {$samples_criteria{$key}{$output_headers[$id_col_array]} eq $_} @excluded) {next;} 
			if ($samples_criteria{$key}{$output_headers[$aff_col_array]} =~ /$type1/i) {
				push @criteria1, $key;
				$count_crit1+=1;
				push @criteria_group, $key;
				$criteria_count+=1;
				push @check_crit_cnv_id, $samples_criteria{$key}{CNV_ID};
				print $OUT_temp2 $type1."\t".$samples_criteria{$key}{$output_headers[$size_col_array]}."\n";
				if (grep {$samples_criteria{$key}{$output_headers[$id_col_array]} eq $_} @ids_crit_total) {
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Size} = $samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Size}+$samples_criteria{$key}{$output_headers[$size_col_array]};
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Genes} = $samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Genes}.",".$samples_criteria{$key}{$output_headers[$gene_col_array]};
						$numcnvs1+=1;
						}
					else {
						push @ids_crit_total, $samples_criteria{$key}{$output_headers[$id_col_array]};
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{ID} = $samples_criteria{$key}{$output_headers[$id_col_array]};
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Pheno} = $samples_criteria{$key}{$output_headers[$aff_col_array]};
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Type} = $samples_criteria{$key}{$output_headers[$type_col_array]};
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Gender} = $samples_criteria{$key}{$output_headers[$gender_col_array]};
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Size} = $samples_criteria{$key}{$output_headers[$size_col_array]};
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Genes} = $samples_criteria{$key}{$output_headers[$gene_col_array]};
						$numcnvs1+=1;
						$totalcritind1+=1;
						}
					}
			if ($samples_criteria{$key}{$output_headers[$aff_col_array]} =~ /$type2/i) {
				push @criteria2, $key; 
				$count_crit2+=1;
				push @criteria_group, $key;
				$criteria_count+=1;
				push @check_crit_cnv_id, $samples_criteria{$key}{CNV_ID};
				print $OUT_temp2 $type2."\t".$samples_criteria{$key}{$output_headers[$size_col_array]}."\n";
				if (grep {$samples_criteria{$key}{$output_headers[$id_col_array]} eq $_} @ids_crit_total) {
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Size} = $samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Size}+$samples_criteria{$key}{$output_headers[$size_col_array]};
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Genes} = $samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Genes}.",".$samples_criteria{$key}{$output_headers[$gene_col_array]};
						$numcnvs2+=1;
						}
					else {
						push @ids_crit_total, $samples_criteria{$key}{$output_headers[$id_col_array]};
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{ID} = $samples_criteria{$key}{$output_headers[$id_col_array]};
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Pheno} = $samples_criteria{$key}{$output_headers[$aff_col_array]};
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Type} = $samples_criteria{$key}{$output_headers[$type_col_array]};
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Gender} = $samples_criteria{$key}{$output_headers[$gender_col_array]};
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Size} = $samples_criteria{$key}{$output_headers[$size_col_array]};
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Genes} = $samples_criteria{$key}{$output_headers[$gene_col_array]};
						$numcnvs2+=1;
						$totalcritind2+=1;
						}
					}
			}
		foreach my $id (@ids_crit_total) {
			print $OUT_total2 "$samples_crit_total{$id}{Pheno}\t$samples_crit_total{$id}{Size}\n";
			}	
		$total_crit_mean_diff = rscript_calc_means($total_crit_outfile,$type1,$type2);
		print $OUTI "The CNV data for the individual totals for all criteria available CNVs.\n";
		($info_data_criteria) = rscript_extract_info($total_crit_outfile,$outfile_final_info);
		print $OUTF "Criteria CNVs\tComparing $type1 and $type2\t$total_crit_mean_diff\t$count_crit1\t$count_crit2\t$totalcritind1\t$totalcritind2\n";
		print $OUT "0\tCriteria CNVs\tComparing $type1 and $type2\t$total_crit_mean_diff\n";
		print $OUTC "0\tCriteria CNVs\tComparing $type1 and $type2\t$total_crit_mean_diff\n";
		print $OUTE "0\tCriteria CNVs\tComparing $type1 and $type2\t$total_crit_mean_diff\n";
		print $OUTEG "0\tCriteria CNVs\tComparing $type1 and $type2\t$total_crit_mean_diff\n";
		print $OUTG "0\tCriteria CNVs\tComparing $type1 and $type2\t$total_crit_mean_diff\n";
		$comparison_text = "Comparing $type1 and $type2";
		if ($number_for_groups) {
			print $OUTN "0\tCriteria CNVs $num_gp1 vs $num_gp2\tComparing $type1 and $type2\t$numbers_mean\n";
			print $OUTF "Criteria CNVs $num_gp1 vs $num_gp2\tComparing $type1 and $type2\t$numbers_mean\n";}
		print $OUTFT "Sample_ID\tPhenotype\tSample_Type\tGender\tCNV_Size\tGenes\n";
		foreach my $ids (@ids_crit_total){
			print $OUTFT "$samples_crit_total{$ids}{ID}\t$samples_crit_total{$ids}{Pheno}\t$samples_crit_total{$ids}{Type}\t$samples_crit_total{$ids}{Gender}\t$samples_crit_total{$ids}{Size}\t$samples_crit_total{$ids}{Genes}\n";
			}
		}
		
	elsif ($type) {
		print "Comparing the types $type.\n";
		print $OUT_temp "Type\tSize\n";
		print $OUT_temp1 "Type\tSize\n";
		print $OUT_temp2 "Type\tSize\n";
		print $OUT_total "Type\tSize\n";
		print $OUT_total1 "Type\tSize\n";
		print $OUT_total2 "Type\tSize\n";
		foreach my $key (@sample_keys){
			if (grep {$samples{$key}{$output_headers[$id_col_array]} eq $_} @excluded) {next;} 
			if ($samples{$key}{$output_headers[$type_col1]} =~ /$type1/i) {
				push @group1, $key;
				$count_group1+=1;
				print $OUT_temp $type1."\t".$samples{$key}{$output_headers[$size_col_array]}."\n";
				if (grep {$samples{$key}{$output_headers[$id_col_array]} eq $_} @ids_total) {
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Size} = $samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Size}+$samples{$key}{$output_headers[$size_col_array]};
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Genes} = $samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Genes}.",".$samples{$key}{$output_headers[$gene_col_array]};
						}
					else {
						push @ids_total, $samples{$key}{$output_headers[$id_col_array]};
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{ID} = $samples{$key}{$output_headers[$id_col_array]};
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Pheno} = $samples{$key}{$output_headers[$aff_col_array]};
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Type} = $samples{$key}{$output_headers[$type_col_array]};
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Gender} = $samples{$key}{$output_headers[$gender_col_array]};
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Size} = $samples{$key}{$output_headers[$size_col_array]};
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Genes} = $samples{$key}{$output_headers[$gene_col_array]};
						$totalind1+=1;
						}
					}
			if ($samples{$key}{$output_headers[$type_col2]} =~ /$type2/i) {
				push @group2, $key; 
				$count_group2+=1;
				print $OUT_temp $type2."\t".$samples{$key}{$output_headers[$size_col_array]}."\n";
				if (grep {$samples{$key}{$output_headers[$id_col_array]} eq $_} @ids_total) {
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Size} = $samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Size}+$samples{$key}{$output_headers[$size_col_array]};
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Genes} = $samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Genes}.",".$samples{$key}{$output_headers[$gene_col_array]};
						}
					else {
						push @ids_total, $samples{$key}{$output_headers[$id_col_array]};
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{ID} = $samples{$key}{$output_headers[$id_col_array]};
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Pheno} = $samples{$key}{$output_headers[$aff_col_array]};
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Type} = $samples{$key}{$output_headers[$type_col_array]};
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Gender} = $samples{$key}{$output_headers[$gender_col_array]};
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Size} = $samples{$key}{$output_headers[$size_col_array]};
						$samples_total{$samples{$key}{$output_headers[$id_col_array]}}{Genes} = $samples{$key}{$output_headers[$gene_col_array]};
						$totalind2+=1;
						}
					}
			}
		foreach my $id (@ids_total) {
			print $OUT_total "$samples_total{$id}{Type}\t$samples_total{$id}{Size}\n";
			}
		$total_mean_diff = rscript_calc_means($total_outfile,$type1,$type2);
		print $OUTI "The CNV data for the individual totals for all available CNVs.\n";
		($info_data_all) = rscript_extract_info($total_outfile,$outfile_final_info);
		#print $LOG1 "The CNV data for the individual totals for all available CNVs.\n$info_data_all\n";
		print $OUTF "All CNVs\tComparing $type1 and $type2\t$total_mean_diff\t$count_group1\t$count_group2\t$totalind1\t$totalind2\n";
		foreach my $key (@exon_keys){
			$samples_exon{$key}{CNV_ID} = $samples_exon{$key}{$output_headers[$id_col_array]}."_".$samples_exon{$key}{$output_headers[$size_col_array]};
			if (grep {$samples_exon{$key}{$output_headers[$id_col_array]} eq $_} @excluded) {next;} 
			if ($samples_exon{$key}{$output_headers[$type_col1]} =~ /$type1/i) {
				push @exon1, $key;
				$count_exon1+=1;
				print $OUT_temp1 $type1."\t".$samples_exon{$key}{$output_headers[$size_col_array]}."\n";
				if (grep {$samples_exon{$key}{$output_headers[$id_col_array]} eq $_} @ids_exon_total) {
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Size} = $samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Size}+$samples_exon{$key}{$output_headers[$size_col_array]};
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Genes} = $samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Genes}.",".$samples_exon{$key}{$output_headers[$gene_col_array]};
						}
					else {
						push @ids_exon_total, $samples_exon{$key}{$output_headers[$id_col_array]};
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{ID} = $samples_exon{$key}{$output_headers[$id_col_array]};
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Pheno} = $samples_exon{$key}{$output_headers[$aff_col_array]};
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Type} = $samples_exon{$key}{$output_headers[$type_col_array]};
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Gender} = $samples_exon{$key}{$output_headers[$gender_col_array]};
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Size} = $samples_exon{$key}{$output_headers[$size_col_array]};
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Genes} = $samples_exon{$key}{$output_headers[$gene_col_array]};
						$totalexonind1+=1;
						}
					}
			if ($samples_exon{$key}{$output_headers[$type_col2]} =~ /$type2/i) {
				push @exon2, $key; 
				$count_exon2+=1;
				print $OUT_temp1 $type2."\t".$samples_exon{$key}{$output_headers[$size_col_array]}."\n";
				if (grep {$samples_exon{$key}{$output_headers[$id_col_array]} eq $_} @ids_exon_total) {
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Size} = $samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Size}+$samples_exon{$key}{$output_headers[$size_col_array]};
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Genes} = $samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Genes}.",".$samples_exon{$key}{$output_headers[$gene_col_array]};
						}
					else {
						push @ids_exon_total, $samples_exon{$key}{$output_headers[$id_col_array]};
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{ID} = $samples_exon{$key}{$output_headers[$id_col_array]};
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Pheno} = $samples_exon{$key}{$output_headers[$aff_col_array]};
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Type} = $samples_exon{$key}{$output_headers[$type_col_array]};
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Gender} = $samples_exon{$key}{$output_headers[$gender_col_array]};
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Size} = $samples_exon{$key}{$output_headers[$size_col_array]};
						$samples_exon_total{$samples_exon{$key}{$output_headers[$id_col_array]}}{Genes} = $samples_exon{$key}{$output_headers[$gene_col_array]};
						$totalexonind2+=1;
						}
					}
			}
		foreach my $id (@ids_exon_total) {
			print $OUT_total1 "$samples_exon_total{$id}{Type}\t$samples_exon_total{$id}{Size}\n";
			}
		$total_exon_mean_diff = rscript_calc_means($total_exon_outfile,$type1,$type2);
		print $OUTI "The CNV data for the individual totals for all exonarray available CNVs.\n";
		($info_data_exon) = rscript_extract_info($total_exon_outfile,$outfile_final_info);
		#print $LOG1 "The CNV data for the individual totals for all exonarray available CNVs.\n$info_data_exon\n";
		print $OUTF "Exonarray CNVs\tComparing $type1 and $type2\t$total_exon_mean_diff\t$count_exon1\t$count_exon2\t$totalexonind1\t$totalexonind2\n";
		foreach my $key (@criteria_keys){
			$samples_criteria{$key}{CNV_ID} = $samples_criteria{$key}{$output_headers[$id_col_array]}."_".$samples_criteria{$key}{$output_headers[$size_col_array]};
			if (grep {$samples_criteria{$key}{$output_headers[$id_col_array]} eq $_} @excluded) {next;} 
			if ($samples_criteria{$key}{$output_headers[$type_col1]} =~ /$type1/i) {
				push @criteria1, $key;
				$count_crit1+=1;
				push @criteria_group, $key;
				$criteria_count+=1;
				push @check_crit_cnv_id, $samples_criteria{$key}{CNV_ID};
				print $OUT_temp2 $type1."\t".$samples_criteria{$key}{$output_headers[$size_col_array]}."\n";
				if (grep {$samples_criteria{$key}{$output_headers[$id_col_array]} eq $_} @ids_crit_total) {
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Size} = $samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Size}+$samples_criteria{$key}{$output_headers[$size_col_array]};
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Genes} = $samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Genes}.",".$samples_criteria{$key}{$output_headers[$gene_col_array]};
						$numcnvs1+=1;
						}
					else {
						push @ids_crit_total, $samples_criteria{$key}{$output_headers[$id_col_array]};
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{ID} = $samples_criteria{$key}{$output_headers[$id_col_array]};
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Pheno} = $samples_criteria{$key}{$output_headers[$aff_col_array]};
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Type} = $samples_criteria{$key}{$output_headers[$type_col_array]};
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Gender} = $samples_criteria{$key}{$output_headers[$gender_col_array]};
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Size} = $samples_criteria{$key}{$output_headers[$size_col_array]};
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Genes} = $samples_criteria{$key}{$output_headers[$gene_col_array]};
						$numcnvs1+=1;
						$totalcritind1+=1;
						}
					}
			if ($samples_criteria{$key}{$output_headers[$type_col2]} =~ /$type2/i) {
				push @criteria2, $key; 
				$count_crit2+=1;
				push @criteria_group, $key;
				$criteria_count+=1;
				push @check_crit_cnv_id, $samples_criteria{$key}{CNV_ID};
				print $OUT_temp2 $type2."\t".$samples_criteria{$key}{$output_headers[$size_col_array]}."\n";
				if (grep {$samples_criteria{$key}{$output_headers[$id_col_array]} eq $_} @ids_crit_total) {
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Size} = $samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Size}+$samples_criteria{$key}{$output_headers[$size_col_array]};
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Genes} = $samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Genes}.",".$samples_criteria{$key}{$output_headers[$gene_col_array]};
						$numcnvs2+=1;
						}
					else {
						push @ids_crit_total, $samples_criteria{$key}{$output_headers[$id_col_array]};
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{ID} = $samples_criteria{$key}{$output_headers[$id_col_array]};
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Pheno} = $samples_criteria{$key}{$output_headers[$aff_col_array]};
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Type} = $samples_criteria{$key}{$output_headers[$type_col_array]};
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Gender} = $samples_criteria{$key}{$output_headers[$gender_col_array]};
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Size} = $samples_criteria{$key}{$output_headers[$size_col_array]};
						$samples_crit_total{$samples_criteria{$key}{$output_headers[$id_col_array]}}{Genes} = $samples_criteria{$key}{$output_headers[$gene_col_array]};
						$numcnvs2+=1;
						$totalcritind2+=1;
						}
					}
			}
		foreach my $id (@ids_crit_total) {
			print $OUT_total2 "$samples_crit_total{$id}{Type}\t$samples_crit_total{$id}{Size}\n";
			}	
		$total_crit_mean_diff = rscript_calc_means($total_crit_outfile,$type1,$type2);
		print $OUTI "The CNV data for the individual totals for all criteria available CNVs.\n";
		($info_data_criteria) = rscript_extract_info($total_crit_outfile,$outfile_final_info);
		print $OUTF "Criteria CNVs\tComparing $type1 and $type2\t$total_crit_mean_diff\t$count_crit1\t$count_crit2\t$totalcritind1\t$totalcritind2\n";
		print $OUT "0\tCriteria CNVs\tComparing $type1 and $type2\t$total_crit_mean_diff\n";
		print $OUTC "0\tCriteria CNVs\tComparing $type1 and $type2\t$total_crit_mean_diff\n";
		print $OUTE "0\tCriteria CNVs\tComparing $type1 and $type2\t$total_crit_mean_diff\n";
		print $OUTEG "0\tCriteria CNVs\tComparing $type1 and $type2\t$total_crit_mean_diff\n";
		print $OUTG "0\tCriteria CNVs\tComparing $type1 and $type2\t$total_crit_mean_diff\n";
		$comparison_text = "Comparing $type1 and $type2";
		if ($number_for_groups) {
			print $OUTN "0\tCriteria CNVs $num_gp1 vs $num_gp2\tComparing $type1 and $type2\t$numbers_mean\n";
			print $OUTF "0\tCriteria CNVs $num_gp1 vs $num_gp2\tComparing $type1 and $type2\t$numbers_mean\n";}
		print $OUTFT "Sample_ID\tPhenotype\tSample_Type\tGender\tCNV_Size\tGenes\n";
		foreach my $ids (@ids_crit_total){
			print $OUTFT "$samples_crit_total{$ids}{ID}\t$samples_crit_total{$ids}{Pheno}\t$samples_crit_total{$ids}{Type}\t$samples_crit_total{$ids}{Gender}\t$samples_crit_total{$ids}{Size}\t$samples_crit_total{$ids}{Genes}\n";
			}
		}
#	}	
### finished getting original p-values and mean difference from files
	
my ($total_all_count,$total_exon_count,$total_criteria_count);
my (@total_all,@total_exon,@total_criteria);
push @total_all, @group1; 
push @total_all, @group2; 
$total_all_count = $count_group1 + $count_group2;
push @total_exon, @exon1; 
push @total_exon, @exon2; 
$total_exon_count = $count_exon1 + $count_exon2;
push @total_criteria, @criteria1; 
push @total_criteria, @criteria2; 
$total_criteria_count = $count_crit1 + $count_crit2;
#print "All the cohort samples keys are:  @total_all.\n";
#print "All the exonarray keys are:  @total_exon.\n";
print "All the criteria keys are:  @total_criteria.\n";	

my $perm_mean_diff;
my ($total_perm_p_value,$total_perm_mean_diff,$final_total_perm_p_value);
my ($count_for_totals,$total_count_perm_p_value,$total_count_perm_mean_diff);
my $final_mean_count_perm_p_value;
my (@total_perm_ids,@total_perm_count_ids);
my %total_perm_samples;
my (%count_for_totals_ids,%total_perm_count_samples);
my ($count_g_perm1)=0;
my $rand_ref1;
my @random_nums1;
my $rand_ref2;
my @random_nums2;


	
#my @sorted_group1 = sort {lc $a cmp lc $b} keys %group1;

my $number_iterations = 1 ;
if ($permutations) {
$starttime = localtime();
print $LOG1 "$starttime:  Starting the permutations for $permutations iterations.\n";
print $OUTF "#Permutations\tComparison\tNum Pick1\tNum Select1\tNum Pick2\tNum Select2\tOriginal Mean Difference\tFinal Total CNV Mean Permutation p-value\n";
	

$number_iterations = $permutations;

print "\n\nStartng permutations for CNV sizes for $number_iterations iterations.\n\n";
	for (my $loop = 1; $loop <= $number_iterations; $loop++){
	
	print "Permutation number:  $loop.\n";
	$starttime = localtime();
	print $LOG1 "$starttime:  Permutation number: $loop\t";
	
	print $OUT $loop."\tAll CNVs available within the groups\t".$comparison_text."\t";
	print $OUTE $loop."\tExonarray CNVs all available within the groups\t".$comparison_text."\t";
	print $OUTEG $loop."\tExonarray CNVs group specific\t".$comparison_text."\t";
	print $OUTG $loop."\tAll CNVs group specific\t".$comparison_text."\t";
	print $OUTC $loop."\tCriteria CNVs randomized\t".$comparison_text."\t";

	my $temp_total_perm = "${output_name}_mean_diff_perm_total.txt";
	open my $OUT_temp_perm_total, '>', $temp_total_perm;
	my $temp_total_perm_count = "${output_name}_mean_diff_perm_total_count.txt";
	open my $OUT_temp_perm_total_count, '>', $temp_total_perm_count;
	print $OUT_temp_perm_total "Type\tSize\n";
	print $OUT_temp_perm_total_count "Type\tSize\n";


	#print "To look at each group and randomly select number for 2 groups to perform tests on.\n";
	
	$rand_ref1 = &random_generator($criteria_count,$total_all_count);
	@random_nums1 = @$rand_ref1;
	$count_for_totals = 0;
	print "Count for totals is:  $count_for_totals.\n";
	for (my $i=0; $i<$criteria_count; $i++){
		if ($i >= $count_crit1) {
			#print "$i\t$random_nums1[$i]\t$total_all[$random_nums1[$i]]\n";
			#print "$samples{$total_all[$random_nums1[$i]]}{$output_headers[$size_col_array]}.\n";
			##print "$samples{$total_all[$random_nums1[$i]]}{$output_headers[$id_col_array]}.\n";
			if ((grep {$samples{$total_all[$random_nums1[$i]]}{$output_headers[$id_col_array]} eq $_} @total_perm_count_ids) &&	($count_for_totals_ids{$samples{$total_all[$random_nums1[$i]]}{$output_headers[$id_col_array]}}{2})) {		
					$total_perm_count_samples{$count_for_totals_ids{$samples{$total_all[$random_nums1[$i]]}{$output_headers[$id_col_array]}}{2}}{Size} = $total_perm_count_samples{$count_for_totals_ids{$samples{$total_all[$random_nums1[$i]]}{$output_headers[$id_col_array]}}{2}}{Size}+$samples{$total_all[$random_nums1[$i]]}{$output_headers[$size_col_array]};
					##print $count_for_totals_ids{$samples{$total_all[$random_nums1[$i]]}{$output_headers[$id_col_array]}}{2}."\n";
				}
			else {
				push @total_perm_count_ids, $samples{$total_all[$random_nums1[$i]]}{$output_headers[$id_col_array]};
				$count_for_totals_ids{$samples{$total_all[$random_nums1[$i]]}{$output_headers[$id_col_array]}}{2} = $count_for_totals;
				$total_perm_count_samples{$count_for_totals}{Size} = $samples{$total_all[$random_nums1[$i]]}{$output_headers[$size_col_array]};
				$total_perm_count_samples{$count_for_totals}{Type} = 2;
				##print "$samples{$total_all[$random_nums1[$i]]}{$output_headers[$id_col_array]}\tType 2\t$count_for_totals.\n";
				$count_for_totals+=1;
				}
			}
		else {
			#print "$i\t$random_nums1[$i]\t$total_all[$random_nums1[$i]]\n";
			#print "$samples{$total_all[$random_nums1[$i]]}{$output_headers[$size_col_array]}.\n";
			#print "$samples{$random_nums1[$i]}{$output_headers[$size_col_array]}.\n";
			if ((grep {$samples{$total_all[$random_nums1[$i]]}{$output_headers[$id_col_array]} eq $_} @total_perm_count_ids) && ($count_for_totals_ids{$samples{$total_all[$random_nums1[$i]]}{$output_headers[$id_col_array]}}{1})) {		
					$total_perm_count_samples{$count_for_totals_ids{$samples{$total_all[$random_nums1[$i]]}{$output_headers[$id_col_array]}}{1}}{Size} = $total_perm_count_samples{$count_for_totals_ids{$samples{$total_all[$random_nums1[$i]]}{$output_headers[$id_col_array]}}{1}}{Size}+$samples{$total_all[$random_nums1[$i]]}{$output_headers[$size_col_array]};
					##print "$samples{$total_all[$random_nums1[$i]]}{$output_headers[$id_col_array]}\t";
					##print "$count_for_totals_ids{$samples{$total_all[$random_nums1[$i]]}{$output_headers[$id_col_array]}}{1}.\n";
				}
			else {
				push @total_perm_count_ids, $samples{$total_all[$random_nums1[$i]]}{$output_headers[$id_col_array]};
				$count_for_totals_ids{$samples{$total_all[$random_nums1[$i]]}{$output_headers[$id_col_array]}}{1} = $count_for_totals;
				$total_perm_count_samples{$count_for_totals}{Size} = $samples{$total_all[$random_nums1[$i]]}{$output_headers[$size_col_array]};
				$total_perm_count_samples{$count_for_totals}{Type} = 1;
				##print "$samples{$total_all[$random_nums1[$i]]}{$output_headers[$id_col_array]}\tType 1\t$count_for_totals.\n";
				$count_for_totals+=1;
				}
			}
		}
	#print "There are $count_for_totals total individuals used.\n";
	for (my $count = 0; $count < $count_for_totals; $count++) {
		print $OUT_temp_perm_total_count $total_perm_count_samples{$count}{Type}."\t".$total_perm_count_samples{$count}{Size}."\n";
		}
	$total_count_perm_mean_diff = rscript_calc_means($temp_total_perm_count,1,2);
	print $OUT $total_count_perm_mean_diff."\n";
	undef @total_perm_count_ids;
	undef %count_for_totals_ids;
	undef %total_perm_count_samples;
	undef @random_nums1;
	close $OUT_temp_perm_total_count;
	open $OUT_temp_perm_total_count, '>', $temp_total_perm_count;
	print $OUT_temp_perm_total_count "Type\tSize\n";
	print $LOG1 "Finished all\t";
	
	$rand_ref1 = &random_generator($count_crit1,$count_group1);
	@random_nums1 = @$rand_ref1;
	$rand_ref2 = &random_generator($count_crit2,$count_group2);
	@random_nums2 = @$rand_ref2;
	$count_g_perm1 = 0;
	for (my $i=0; $i<$criteria_count; $i++){
		if ($i >= $count_crit1) {
			if (grep {$samples{$group2[$random_nums2[$count_g_perm1]]}{$output_headers[$id_col_array]} eq $_} @total_perm_ids) {
				$total_perm_samples{$samples{$group2[$random_nums2[$count_g_perm1]]}{$output_headers[$id_col_array]}}{Size} = $total_perm_samples{$samples{$group2[$random_nums2[$count_g_perm1]]}{$output_headers[$id_col_array]}}{Size}+$samples{$group2[$random_nums2[$count_g_perm1]]}{$output_headers[$size_col_array]};
				}
			else {
				push @total_perm_ids, $samples{$group2[$random_nums2[$count_g_perm1]]}{$output_headers[$id_col_array]};
				$total_perm_samples{$samples{$group2[$random_nums2[$count_g_perm1]]}{$output_headers[$id_col_array]}}{ID} = $samples{$group2[$random_nums2[$count_g_perm1]]}{$output_headers[$id_col_array]};
				$total_perm_samples{$samples{$group2[$random_nums2[$count_g_perm1]]}{$output_headers[$id_col_array]}}{Size} = $samples{$group2[$random_nums2[$count_g_perm1]]}{$output_headers[$size_col_array]};
				$total_perm_samples{$samples{$group2[$random_nums2[$count_g_perm1]]}{$output_headers[$id_col_array]}}{Type} = 2;
				}
			$count_g_perm1+=1;
			}
		else {
			if (grep {$samples{$group1[$random_nums1[$i]]}{$output_headers[$id_col_array]} eq $_} @total_perm_ids) {
				$total_perm_samples{$samples{$group1[$random_nums1[$i]]}{$output_headers[$id_col_array]}}{Size} = $total_perm_samples{$samples{$group1[$random_nums1[$i]]}{$output_headers[$id_col_array]}}{Size}+$samples{$group1[$random_nums1[$i]]}{$output_headers[$size_col_array]};
				}
			else {
				push @total_perm_ids, $samples{$group1[$random_nums1[$i]]}{$output_headers[$id_col_array]};
				$total_perm_samples{$samples{$group1[$random_nums1[$i]]}{$output_headers[$id_col_array]}}{ID} = $samples{$group1[$random_nums1[$i]]}{$output_headers[$id_col_array]};
				$total_perm_samples{$samples{$group1[$random_nums1[$i]]}{$output_headers[$id_col_array]}}{Size} = $samples{$group1[$random_nums1[$i]]}{$output_headers[$size_col_array]};
				$total_perm_samples{$samples{$group1[$random_nums1[$i]]}{$output_headers[$id_col_array]}}{Type} = 1;
				}
			}
		}
	foreach my $id (@total_perm_ids) {
		print $OUT_temp_perm_total $total_perm_samples{$id}{Type}."\t".$total_perm_samples{$id}{Size}."\n";
		}
	$total_perm_mean_diff = rscript_calc_means($temp_total_perm,1,2);
	print $OUTG $total_perm_mean_diff."\n";
	undef @total_perm_ids;
	undef %total_perm_samples;
	undef @random_nums1;
	undef @random_nums2;
	close $OUT_temp_perm_total;
	open $OUT_temp_perm_total, '>', $temp_total_perm;
	print $OUT_temp_perm_total "Type\tSize\n";
	print $LOG1 "Finished group\t";
	
	my @temp_exon = @total_exon;
	$count_for_totals = 0;
	for (my $i=0; $i<$criteria_count; $i++){
		my $temp_total_exon = scalar @temp_exon;
		my $random_nums1 = int(rand($temp_total_exon));
		print $temp_exon[$random_nums1]." ";
		if ($i >= $count_crit1) {
			if ((grep {$samples_exon{$temp_exon[$random_nums1]}{$output_headers[$id_col_array]} eq $_} @total_perm_count_ids) && ($count_for_totals_ids{$samples_exon{$temp_exon[$random_nums1]}{$output_headers[$id_col_array]}}{2})) {		
					$total_perm_count_samples{$count_for_totals_ids{$samples_exon{$temp_exon[$random_nums1]}{$output_headers[$id_col_array]}}{2}}{Size} = $total_perm_count_samples{$count_for_totals_ids{$samples_exon{$temp_exon[$random_nums1]}{$output_headers[$id_col_array]}}{2}}{Size}+$samples_exon{$temp_exon[$random_nums1]}{$output_headers[$size_col_array]};
				}
			else {
				push @total_perm_count_ids, $samples_exon{$temp_exon[$random_nums1]}{$output_headers[$id_col_array]};
				$count_for_totals_ids{$samples_exon{$temp_exon[$random_nums1]}{$output_headers[$id_col_array]}}{2} = $count_for_totals;
				$total_perm_count_samples{$count_for_totals}{Size} = $samples_exon{$temp_exon[$random_nums1]}{$output_headers[$size_col_array]};
				$total_perm_count_samples{$count_for_totals}{Type} = 2;
				$count_for_totals+=1;
				}
			}
		else {
			if  ((grep {$samples_exon{$temp_exon[$random_nums1]}{$output_headers[$id_col_array]} eq $_} @total_perm_count_ids) && ($count_for_totals_ids{$samples_exon{$temp_exon[$random_nums1]}{$output_headers[$id_col_array]}}{1})) {		
					$total_perm_count_samples{$count_for_totals_ids{$samples_exon{$temp_exon[$random_nums1]}{$output_headers[$id_col_array]}}{1}}{Size} = $total_perm_count_samples{$count_for_totals_ids{$samples_exon{$temp_exon[$random_nums1]}{$output_headers[$id_col_array]}}{1}}{Size}+$samples_exon{$temp_exon[$random_nums1]}{$output_headers[$size_col_array]};
				}
			else {
				push @total_perm_count_ids, $samples_exon{$temp_exon[$random_nums1]}{$output_headers[$id_col_array]};
				$count_for_totals_ids{$samples_exon{$temp_exon[$random_nums1]}{$output_headers[$id_col_array]}}{1} = $count_for_totals;
				$total_perm_count_samples{$count_for_totals}{Size} = $samples_exon{$temp_exon[$random_nums1]}{$output_headers[$size_col_array]};
				$total_perm_count_samples{$count_for_totals}{Type} = 1;
				$count_for_totals+=1;
				}
			}
		splice @temp_exon, $random_nums1, 1;
		}
	print "\n";
	for (my $count=0;$count<$count_for_totals;$count++) {
		print $OUT_temp_perm_total_count $total_perm_count_samples{$count}{Type}."\t".$total_perm_count_samples{$count}{Size}."\n";
		}
	$total_count_perm_mean_diff = rscript_calc_means($temp_total_perm_count,1,2);
	print $OUTE $total_count_perm_mean_diff."\n";
	undef @total_perm_count_ids;
	undef %count_for_totals_ids;
	undef %total_perm_count_samples;
	undef @random_nums1;
	close $OUT_temp_perm_total_count;
	open $OUT_temp_perm_total_count, '>', $temp_total_perm_count;
	print $OUT_temp_perm_total_count "Type\tSize\n";
	print $LOG1 "Finished exonarray\t";
	
	$rand_ref1 = &random_generator($count_crit1,$count_exon1);
	@random_nums1 = @$rand_ref1;
	$rand_ref2 = &random_generator($count_crit2,$count_exon2);
	@random_nums2 = @$rand_ref2;
	$count_g_perm1 = 0;
	for (my $i=0; $i<$criteria_count; $i++){
		if ($i >= $count_crit1) {
			if (grep {$samples_exon{$exon2[$random_nums2[$count_g_perm1]]}{$output_headers[$id_col_array]} eq $_} @total_perm_ids) {
				$total_perm_samples{$samples_exon{$exon2[$random_nums2[$count_g_perm1]]}{$output_headers[$id_col_array]}}{Size} = $total_perm_samples{$samples_exon{$exon2[$random_nums2[$count_g_perm1]]}{$output_headers[$id_col_array]}}{Size} +$samples_exon{$exon2[$random_nums2[$count_g_perm1]]}{$output_headers[$size_col_array]};
				}
			else {
				push @total_perm_ids, $samples_exon{$exon2[$random_nums2[$count_g_perm1]]}{$output_headers[$id_col_array]};
				$total_perm_samples{$samples_exon{$exon2[$random_nums2[$count_g_perm1]]}{$output_headers[$id_col_array]}}{ID} = $samples_exon{$exon2[$random_nums2[$count_g_perm1]]}{$output_headers[$id_col_array]};
				$total_perm_samples{$samples_exon{$exon2[$random_nums2[$count_g_perm1]]}{$output_headers[$id_col_array]}}{Size} = $samples_exon{$exon2[$random_nums2[$count_g_perm1]]}{$output_headers[$size_col_array]};
				$total_perm_samples{$samples_exon{$exon2[$random_nums2[$count_g_perm1]]}{$output_headers[$id_col_array]}}{Type} = 2;
				}
			$count_g_perm1+=1;
			}
		else {
			if (grep {$samples_exon{$exon1[$random_nums1[$i]]}{$output_headers[$id_col_array]} eq $_} @total_perm_ids) {
				$total_perm_samples{$samples_exon{$exon1[$random_nums1[$i]]}{$output_headers[$id_col_array]}}{Size} = $total_perm_samples{$samples_exon{$exon1[$random_nums1[$i]]}{$output_headers[$id_col_array]}}{Size} + $samples_exon{$exon1[$random_nums1[$i]]}{$output_headers[$size_col_array]};
			}
			else {
				push @total_perm_ids, $samples_exon{$exon1[$random_nums1[$i]]}{$output_headers[$id_col_array]};
				$total_perm_samples{$samples_exon{$exon1[$random_nums1[$i]]}{$output_headers[$id_col_array]}}{ID} = $samples_exon{$exon1[$random_nums1[$i]]}{$output_headers[$id_col_array]};
				$total_perm_samples{$samples_exon{$exon1[$random_nums1[$i]]}{$output_headers[$id_col_array]}}{Size} = $samples_exon{$exon1[$random_nums1[$i]]}{$output_headers[$size_col_array]};
				$total_perm_samples{$samples_exon{$exon1[$random_nums1[$i]]}{$output_headers[$id_col_array]}}{Type} = 1;
				}
			}
		}
	foreach my $id (@total_perm_ids) {
		print $OUT_temp_perm_total $total_perm_samples{$id}{Type}."\t".$total_perm_samples{$id}{Size}."\n";
		}
	$total_perm_mean_diff = rscript_calc_means($temp_total_perm,1,2);
	print $OUTEG $total_perm_mean_diff."\n";
	undef @total_perm_ids;
	undef %total_perm_samples;
	undef @random_nums1;
	undef @random_nums2;
	close $OUT_temp_perm_total;
	open $OUT_temp_perm_total, '>', $temp_total_perm;
	print $OUT_temp_perm_total "Type\tSize\n";
	print $LOG1 "Finished exon_group\t";
	
	@random_nums1 = shuffle(@criteria_group);  ##shuffling the criteria CNVs
	#print "@random_nums1\n";
	$count_for_totals = 0;
	for (my $i=0; $i<$criteria_count; $i++){
		if ($i >= $count_crit1) {
			#print "$i\t$random_nums1[$i]\t$criteria_group[$random_nums1[$i]]\n";
			#print "$samples_criteria{$criteria_group[$random_nums1[$i]]}{$output_headers[$size_col_array]}.\n";
			##print "$samples_criteria{$random_nums1[$i]}{$output_headers[$size_col_array]}.\n";
			#print $OUT_temp_perm_A "2\t$samples_criteria{$criteria_group[$random_nums1[$i]]}{$output_headers[$size_col_array]}\n";}
			if ((grep {$samples_criteria{$random_nums1[$i]}{$output_headers[$id_col_array]} eq $_} @total_perm_count_ids) && ($count_for_totals_ids{$samples_criteria{$random_nums1[$i]}{$output_headers[$id_col_array]}}{2})) {		
					$total_perm_count_samples{$count_for_totals_ids{$samples_criteria{$random_nums1[$i]}{$output_headers[$id_col_array]}}{2}}{Size} = $total_perm_count_samples{$count_for_totals_ids{$samples_criteria{$random_nums1[$i]}{$output_headers[$id_col_array]}}{2}}{Size}+$samples_criteria{$random_nums1[$i]}{$output_headers[$size_col_array]};
				}
			else {
				push @total_perm_count_ids, $samples_criteria{$random_nums1[$i]}{$output_headers[$id_col_array]};
				$count_for_totals_ids{$samples_criteria{$random_nums1[$i]}{$output_headers[$id_col_array]}}{2} = $count_for_totals;
				$total_perm_count_samples{$count_for_totals}{Size} = $samples_criteria{$random_nums1[$i]}{$output_headers[$size_col_array]};
				$total_perm_count_samples{$count_for_totals}{Type} = 2;
				$count_for_totals+=1;
				}
			}
		else {
			#print "$i\t$random_nums1[$i]\t$criteria_group[$random_nums1[$i]]\n";
			#print "$samples_criteria{$criteria_group[$random_nums1[$i]]}{$output_headers[$size_col_array]}.\n";
			##print "$samples_criteria{$random_nums1[$i]}{$output_headers[$size_col_array]}.\n";
			#print $OUT_temp_perm_A "1\t$samples_criteria{$criteria_group[$random_nums1[$i]]}{$output_headers[$size_col_array]}\n";}
			if ((grep {$samples_criteria{$random_nums1[$i]}{$output_headers[$id_col_array]} eq $_} @total_perm_count_ids) && ($count_for_totals_ids{$samples_criteria{$random_nums1[$i]}{$output_headers[$id_col_array]}}{1})) {		
					$total_perm_count_samples{$count_for_totals_ids{$samples_criteria{$random_nums1[$i]}{$output_headers[$id_col_array]}}{1}}{Size} = $total_perm_count_samples{$count_for_totals_ids{$samples_criteria{$random_nums1[$i]}{$output_headers[$id_col_array]}}{1}}{Size}+$samples_criteria{$random_nums1[$i]}{$output_headers[$size_col_array]};
				}
			else {
				push @total_perm_count_ids, $samples_criteria{$random_nums1[$i]}{$output_headers[$id_col_array]};
				$count_for_totals_ids{$samples_criteria{$random_nums1[$i]}{$output_headers[$id_col_array]}}{1} = $count_for_totals;
				$total_perm_count_samples{$count_for_totals}{Size} = $samples_criteria{$random_nums1[$i]}{$output_headers[$size_col_array]};
				$total_perm_count_samples{$count_for_totals}{Type} = 1;
				$count_for_totals+=1;
				}
			}
		}
	for (my $count=0;$count<$count_for_totals;$count++) {
		print $OUT_temp_perm_total_count $total_perm_count_samples{$count}{Type}."\t".$total_perm_count_samples{$count}{Size}."\n";
		}
	$total_count_perm_mean_diff = rscript_calc_means($temp_total_perm_count,1,2);
	print $OUTC $total_count_perm_mean_diff."\n";
	undef @total_perm_count_ids;
	undef %count_for_totals_ids;
	undef %total_perm_count_samples;
	undef @random_nums1;
	close $OUT_temp_perm_total;
	close $OUT_temp_perm_total_count;
	print $LOG1 "Finished criteria\t";
	
	if ($number_for_groups) {
		print $OUTN "$loop\tCriteria CNVs $num_gp1 vs $num_gp2\t$comparison_text\t";
		open $OUT_temp_perm_total, '>', $temp_total_perm;
		print $OUT_temp_perm_total "Type\tSize\n";
		my $loop_count = $num_gp1 + $num_gp2;
		$rand_ref1 = &random_generator($num_gp1,$count_crit1);
		@random_nums1 = @$rand_ref1;
		$rand_ref2 = &random_generator($num_gp2,$count_crit2);
		@random_nums2 = @$rand_ref2;
		$count_g_perm1=0;
		for (my $i=0; $i<$loop_count; $i++){
			if ($i >= $num_gp1) {
				if (grep {$samples_criteria{$criteria2[$random_nums2[$count_g_perm1]]}{$output_headers[$id_col_array]} eq $_} @total_perm_ids) {
					$total_perm_samples{$samples_criteria{$criteria2[$random_nums2[$count_g_perm1]]}{$output_headers[$id_col_array]}}{Size} = $total_perm_samples{$samples_criteria{$criteria2[$random_nums2[$count_g_perm1]]}{$output_headers[$id_col_array]}}{Size} +$samples_criteria{$criteria2[$random_nums2[$count_g_perm1]]}{$output_headers[$size_col_array]};
					}
				else {
					push @total_perm_ids, $samples_criteria{$criteria2[$random_nums2[$count_g_perm1]]}{$output_headers[$id_col_array]};
					$total_perm_samples{$samples_criteria{$criteria2[$random_nums2[$count_g_perm1]]}{$output_headers[$id_col_array]}}{ID} = $samples_criteria{$criteria2[$random_nums2[$count_g_perm1]]}{$output_headers[$id_col_array]};
					$total_perm_samples{$samples_criteria{$criteria2[$random_nums2[$count_g_perm1]]}{$output_headers[$id_col_array]}}{Size} = $samples_criteria{$criteria2[$random_nums2[$count_g_perm1]]}{$output_headers[$size_col_array]};
					$total_perm_samples{$samples_criteria{$criteria2[$random_nums2[$count_g_perm1]]}{$output_headers[$id_col_array]}}{Type} = 2;
					}
				$count_g_perm1+=1;
				}
			else {
				if (grep {$samples_criteria{$criteria1[$random_nums1[$i]]}{$output_headers[$id_col_array]} eq $_} @total_perm_ids) {
					$total_perm_samples{$samples_criteria{$criteria1[$random_nums1[$i]]}{$output_headers[$id_col_array]}}{Size} = $total_perm_samples{$samples_criteria{$criteria1[$random_nums1[$i]]}{$output_headers[$id_col_array]}}{Size} + $samples_criteria{$criteria1[$random_nums1[$i]]}{$output_headers[$size_col_array]};
					}
				else {
					push @total_perm_ids, $samples_criteria{$criteria1[$random_nums1[$i]]}{$output_headers[$id_col_array]};
					$total_perm_samples{$samples_criteria{$criteria1[$random_nums1[$i]]}{$output_headers[$id_col_array]}}{ID} = $samples_criteria{$criteria1[$random_nums1[$i]]}{$output_headers[$id_col_array]};
					$total_perm_samples{$samples_criteria{$criteria1[$random_nums1[$i]]}{$output_headers[$id_col_array]}}{Size} = $samples_criteria{$criteria1[$random_nums1[$i]]}{$output_headers[$size_col_array]};
					$total_perm_samples{$samples_criteria{$criteria1[$random_nums1[$i]]}{$output_headers[$id_col_array]}}{Type} = 1;
					}
				}
			}
		foreach my $id (@total_perm_ids) {
			print $OUT_temp_perm_total $total_perm_samples{$id}{Type}."\t".$total_perm_samples{$id}{Size}."\n";
			}
		$total_perm_mean_diff = rscript_calc_means($temp_total_perm,1,2);
		print $OUTN $total_perm_mean_diff."\n";
		undef @total_perm_ids;
		undef %total_perm_samples;
		undef @random_nums1;
		undef @random_nums2;
		close $OUT_temp_perm_total;
		print $LOG1 "Finished numbers\t";
		}
		
		qx {rm $temp_total_perm};
		qx {rm $temp_total_perm_count};
		
		print $LOG1 "Finished runs\n";
		
	}
	
	
	if ($total_crit_mean_diff > 0) {
		$final_mean_count_perm_p_value = qx {awk -F'\t' 'NR>2 && \$4>=$total_crit_mean_diff {print \$0}' $outfile_all| wc -l};
		chomp $final_mean_count_perm_p_value;
		}
	else {
		$final_mean_count_perm_p_value = qx {awk -F'\t' 'NR>2 && \$4<=$total_crit_mean_diff {print \$0}' $outfile_all| wc -l};
		chomp $final_mean_count_perm_p_value;
		}
	
	if ($final_mean_count_perm_p_value == 0) {
		print $OUTF "$number_iterations\tAll CNVs available within the groups\t$count_crit1\t$total_all_count\t$count_crit2\t$total_all_count\t$total_crit_mean_diff\t$final_mean_count_perm_p_value\n";
		}
	else {
		$final_mean_count_perm_p_value = $final_mean_count_perm_p_value / $permutations;
		print $OUTF "$number_iterations\tAll CNVs available within the groups\t$count_crit1\t$total_all_count\t$count_crit2\t$total_all_count\t$total_crit_mean_diff\t$final_mean_count_perm_p_value\n";
		}
		
	if ($total_crit_mean_diff > 0) {
		$final_mean_count_perm_p_value = qx {awk -F'\t' 'NR>2 && \$4>=$total_crit_mean_diff {print \$0}' $outfile_group | wc -l};
		chomp $final_mean_count_perm_p_value;
		}
	else {
		$final_mean_count_perm_p_value = qx {awk -F'\t' 'NR>2 && \$4<=$total_crit_mean_diff {print \$0}' $outfile_group | wc -l};
		chomp $final_mean_count_perm_p_value;
		}
	
	if ($final_mean_count_perm_p_value == 0) {
		print $OUTF "$number_iterations\tAll CNVs group specific\t$count_crit1\t$count_group1\t$count_crit2\t$count_group2\t$total_crit_mean_diff\t$final_mean_count_perm_p_value\n";
		}
	else {
		$final_mean_count_perm_p_value = $final_mean_count_perm_p_value / $permutations;
		print $OUTF "$number_iterations\tAll CNVs group specific\t$count_crit1\t$count_group1\t$count_crit2\t$count_group2\t$total_crit_mean_diff\t$final_mean_count_perm_p_value\n";
		}
		
	if ($total_crit_mean_diff > 0) {
		$final_mean_count_perm_p_value = qx {awk -F'\t' 'NR>2 && \$4>=$total_crit_mean_diff {print \$0}' $outfile_exon | wc -l};
		chomp $final_mean_count_perm_p_value;
		}
	else {
		$final_mean_count_perm_p_value = qx {awk -F'\t' 'NR>2 && \$4<=$total_crit_mean_diff {print \$0}' $outfile_exon | wc -l};
		chomp $final_mean_count_perm_p_value;
		}

	if ($final_mean_count_perm_p_value == 0) {
		print $OUTF "$number_iterations\tExonarray CNVs all available within the groups\t$count_crit1\t$total_exon_count\t$count_crit2\t$total_exon_count\t$total_crit_mean_diff\t$final_mean_count_perm_p_value\n";
		}
	else {
		$final_mean_count_perm_p_value = $final_mean_count_perm_p_value / $permutations;
		print $OUTF "$number_iterations\tExonarray CNVs all available within the groups\t$count_crit1\t$total_exon_count\t$count_crit2\t$total_exon_count\t$total_crit_mean_diff\t$final_mean_count_perm_p_value\n";
		}
		
	if ($total_crit_mean_diff > 0) {
		$final_mean_count_perm_p_value = qx {awk -F'\t' 'NR>2 && \$4>=$total_crit_mean_diff {print \$0}' $outfile_exon_group | wc -l};
		chomp $final_mean_count_perm_p_value;
		}
	else {
		$final_mean_count_perm_p_value = qx {awk -F'\t' 'NR>2 && \$4<=$total_crit_mean_diff {print \$0}' $outfile_exon_group | wc -l};
		chomp $final_mean_count_perm_p_value;
		}
	
	if ( $final_mean_count_perm_p_value == 0) {
		print $OUTF "$number_iterations\tExonarray CNVs group specific\t$count_crit1\t$count_exon1\t$count_crit2\t$count_exon2\t$total_crit_mean_diff\t$final_mean_count_perm_p_value\n";
		}
	else {
		$final_mean_count_perm_p_value = $final_mean_count_perm_p_value / $permutations;
		print $OUTF "$number_iterations\tExonarray CNVs group specific\t$count_crit1\t$count_exon1\t$count_crit2\t$count_exon2\t$total_crit_mean_diff\t$final_mean_count_perm_p_value\n";
		}
		
	if ($total_crit_mean_diff > 0) {
		$final_mean_count_perm_p_value = qx {awk -F'\t' 'NR>2 && \$4>=$total_crit_mean_diff {print \$0}' $outCriteria | wc -l};
		chomp $final_mean_count_perm_p_value;
		}
	else {
		$final_mean_count_perm_p_value = qx {awk -F'\t' 'NR>2 && \$4<=$total_crit_mean_diff {print \$0}' $outCriteria | wc -l};
		chomp $final_mean_count_perm_p_value;
		}
	
	if ($final_mean_count_perm_p_value == 0) {
		print $OUTF "$number_iterations\tCriteria CNVs shuffled\t$count_crit1\t$criteria_count\t$count_crit2\t$criteria_count\t$total_crit_mean_diff\t$final_mean_count_perm_p_value\n";
		}
	else {
		$final_mean_count_perm_p_value = $final_mean_count_perm_p_value / $permutations;
		print $OUTF "$number_iterations\tCriteria CNVs shuffled\t$count_crit1\t$criteria_count\t$count_crit2\t$criteria_count\t$total_crit_mean_diff\t$final_mean_count_perm_p_value\n";
		#print "Looking at mean $total_crit_mean_diff, the final mean p is $final_mean_count_perm_p_value.\n";
		}
	
	if ($number_for_groups) {
		$final_mean_count_perm_p_value = qx {awk -F'\t' 'NR>2 && \$4>=$numbers_mean {print \$0}' $outfile_numbers | wc -l};
		chomp $final_mean_count_perm_p_value;
		if ($final_mean_count_perm_p_value == 0) {
			print $OUTF "$number_iterations\tCriteria CNVs $num_gp1 vs $num_gp2\t$num_gp1\t$count_crit1\t$num_gp2\t$count_crit2\t$numbers_mean\t$final_mean_count_perm_p_value\n";
			}
		else {
			$final_mean_count_perm_p_value = $final_mean_count_perm_p_value / $permutations;
			print $OUTF "$number_iterations\tCriteria CNVs $num_gp1 vs $num_gp2\t$num_gp1\t$count_crit1\t$num_gp2\t$count_crit2\t$numbers_mean\t$final_mean_count_perm_p_value\n";
			}
		}
	
	print $LOG1 "Finished Final p-value calculations for $permutations iterations run on the different potential mean differences\n\n";
	print "Finished Final p-value calculations for $permutations iterations run on the different potential mean differences\n\n";
	
		
}
else {
	##if not doing permutations, close and remove permutation files
	close $OUT; 
	close $OUTC; 
	close $OUTG; 
	close $OUTE; 
	close $OUTEG; 
	qx {rm $outfile_all};
	qx {rm $outfile_exon};
	qx {rm $outfile_group};
	qx {rm $outfile_exon_group};
	qx {rm $outCriteria};
	
	if ($number_for_groups) {
		close $OUTN; 
		qx {rm $outfile_numbers};
		}
	}


##close temp and total outfiles
close $OUT_temp; close $OUT_temp1; close $OUT_temp2; close $OUT_total; close $OUT_total1; close $OUT_total2; 

#close outfiles
if ($permutations) {
	if ($number_for_groups) {
		close $OUTN; 
		}
	close $OUT; 
	close $OUTE; 
	close $OUTEG; 
	close $OUTG;
	close $OUTC; 
	}	 
close $OUTF; 
				

	$starttime = localtime();
	print $LOG1 "$starttime:  Finished running program $0.\n";
	my ($minutes,$seconds) = runtime($original_starttime);
	print $LOG1 "\nProgramme took $minutes min $seconds s to run\n";

#close INFILE;
close $LOG1;



#################################################################################
##  SUBROUTINES ##
#################################################################################


#random generator to use for iterations
sub random_generator {	
  my ($howmany, $range) = @_; 	#number_to_pick, $total_to_pick_from
  my @num_array = ();
  my $random_number;
  print "You want ",$howmany," numbers to be generated from 0 to ",$range,".\n\n";
  
  for (my $i=0; $i<$howmany; $i++){
	#print "Picking number for #$i.\n";
	$random_number = int(rand($range));
	if (grep (/$random_number/,@num_array)) {
		#print "$random_number already seen.\n";
		#print "Need to pick a new number for #$i:  Going back to ";
		$i--;
		#print "#$i before the loop.\n";
		next;
		}
	else {
		#print "This random number is $random_number.\n";
		push @num_array, $random_number;
		}
	}
#  print "Here are the random numbers:  \n\n";
  undef $howmany;
  undef $range;
  undef $random_number;
 # print "@num_array\n";
  return (\@num_array);
  undef @num_array;
}


#running R script
sub rscript_calc_means {
	my ($outfile,$first,$second) = @_;
	my @means;
	my $first_group;
	my $second_group;
	print "Running the Rscript to calculate the means on file $outfile.\n";
	my $line = qx {Rscript $rscript_file $outfile $first $second 2>/dev/null | grep -i "mean"};
	chomp $line;
	@means = split ('\t',$line);
	@means = map {$_ =~ s/^\s+//; $_} @means;
	@means = map {$_ =~ s/\s+$//; $_} @means;
	$first_group = $means[1];
	$second_group = $means[2];
	print "The mean for $first is $first_group.\n";
	print "The mean for $second is $second_group.\n";
	my $mean_diff = $first_group - $second_group;
	print "The mean difference is $mean_diff.\n";
	return ($mean_diff);
}

sub rscript_extract_info {
	my ($outfile, $file) = @_;
	#qx {Rscript $rscript_extract_info $outfile $file 2>/dev/null};
	my $info = qx {Rscript $rscript_extract_info $outfile $file 2>/dev/null | sed -n '/N/,\$p'};
	chomp $info;
	return ($info);
}


#Calculates the length of time a process has been running for
sub runtime {
	my $starttime = $_[0];
	my $endtime = time;
	my $lengthtime = $endtime - $starttime;
	my $minutes = floor($lengthtime/60);
	my $seconds = sprintf("%02d", ($lengthtime - $minutes*60));
	#unless ($_[1] eq "noprint"){
	#	print "\nProgramme took $minutes min $seconds s to run\n";
	#}
	return ($minutes, $seconds);
}



	#================================================================================================
	### DESCRIPTION: To make a file with the median scores of expression data for each brain region/time for a list of genes 
	### Please provide the list of genes, the expression data file, the column data sheet, the genecode file, and output file name
	### EX: <genes_of_interest><expression_data_matrix><columns_data4.txt><genecode.txt><output_name>
	### OUTPUT: Will provide a new text file containing the median expression values for the list of genes using regions of choice
	### data will be separated by hemisphere, for each time period 3-15, for each brain region 
	#==================================================================================================\n";
sub pulling_out_exonarray_data {
	my ($gene_list, $expression_data, $columns_data, $gene_code, $output_name) = @_;
	
	my $outfile_pull_data;
	my $OUT_pd;
	my $outfile_list;
	my $OUTL2;
	my @genes;
	my @gene_codes_list;
	my @gene_names_list;
	my $final_gene_count;
	my $gene_count=0;
	my @columns_num;
	my $column;
	my $column_num;
	my @brain_regions;
	my @time_periods;
	my @all_genes;
	my $temp_list;
	my $temp_outfile;
	my @hemis = ('L','R');
	my $hemi_column = qx {head -1 $columns_data | awk '{for (i=1;i<NF;i++) if(\$i~/Hemisphere/){print i}}'};
	chomp $hemi_column;
	my $time_column = qx {head -1 $columns_data | awk '{for (i=1;i<NF;i++) if(\$i~/Time/){print i}}'};
	chomp $time_column;
	my $ref_names;
	my $ref_codes;

	$temp_list = qx {awk '{print \$0}' $gene_list};
	chomp $temp_list;
	@genes = split('\n',$temp_list);
	$starttime = localtime();
	print "This list has ",scalar(@genes)," number of genes to look up in it.\n";
	print $LOG1 "$starttime: This list $gene_list has ",scalar(@genes)," number of genes to look up in it.\n";

	(my $stem) = $gene_list =~ /(.+).txt/;
	$temp_outfile = "${stem}_temp.txt";
	print "Created temporary file for exonarray values of genes in $temp_outfile\n";

	my $temp_code;
	foreach my $list_genes (@genes) {
		chomp $list_genes;
#		print "**$list_genes**\n";
		$temp_code = qx {awk '\$3=="$list_genes" {print \$1}' $gene_code};
		chomp $temp_code;
#		print "**$temp_code**\n";
		if ($temp_code eq "" || $temp_code eq "\n") {next;}
		elsif ($temp_code =~ /\s/) {next;}
		elsif (grep (/$temp_code/,@gene_codes_list)) {next;}
		else {push @gene_codes_list, $temp_code;}
		}
	#print "@gene_codes_list\n";
	$starttime = localtime();
	print "There were ",scalar(@gene_codes_list)," genes codes found from the original list to put into the temporary file.\n\n";
	print $LOG1 "$starttime: There were ",scalar(@gene_codes_list)," genes codes found from the original list to put into the temporary file.\n\n";

	$outfile_list = "${output_name}_final_genelist.txt"; 
	open $OUTL2, '>', $outfile_list;
	print $LOG1 "Created $outfile_list to contain the final genelist of the exonarray data for the genes available.\n";
	my $temp1;
#	print "My gene codes file is ",$gene_code,".\n\n";
	print "Finding genes of interests names!\n\n";
	foreach my $temp_genes (@gene_codes_list) {
		#print "$temp_genes\n";
		$temp1 = qx {awk '\$1=="$temp_genes" {print \$3}' $gene_code};
		chomp $temp1;
		#print "$temp1\n";
		push @gene_names_list, $temp1;
	}
#	@gene_names_list = @$ref_names;
	foreach my $name (@gene_names_list) {
		print $OUTL2 $name."\n";
		}
	close $OUTL2;

	foreach my $code (@gene_codes_list){
		if ($code == $gene_codes_list[0]){
			qx {awk '\$1==$code {print \$0}' $expression_data > $temp_outfile};}
		else {qx {awk '\$1==$code {print \$0}' $expression_data >> $temp_outfile};}
		}

	#open temp file to look at each row
	open INFILE , '<' , $temp_outfile or die "Could not open file $temp_outfile\n$!";

	#determine what type of regions to use
	print "Looking at all brain regions.\n\n";
	@brain_regions = ('OFC','DFC','VFC','MFC','M1C','S1C','IPC','A1C','STC','ITC','V1C','HIP','AMY','STR','MD','CBC');
	@time_periods = (3..15);
	$column = qx {head -1 $columns_data | awk '{for (i=1;i<NF;i++) if(\$i~/Individual/){print i}}'};
	chomp $column;
	$outfile_pull_data = "${output_name}.txt";
	open $OUT_pd, '>', $outfile_pull_data;
	print $LOG1 "Created $outfile_pull_data to contain the exonarray data for the genes available.\n\n";

	###Start making file
	print $OUT_pd "Gene\t";
	foreach my $time (@time_periods){
		foreach my $region (@brain_regions){
			foreach my $hemis (@hemis){
				print $OUT_pd $time."_".$region."_".$hemis."\t";
				}
			}
		}
		print $OUT_pd "\n";

	print "The columns will be found in $columns_data file\n";
	while (my $line = <INFILE>){
		$line =~ s/[\n\r]+//g;
		my @tab = split /\t/,$line;
		print $OUT_pd $gene_names_list[$gene_count];
		foreach my $time (@time_periods){
			print "This time period is $time\n";
			foreach my $brain_region (@brain_regions){
				foreach my $hemis (@hemis){
					$column_num = qx {awk '\$${time_column}=="$time" && \$${hemi_column}=="$hemis" && \$${column}=="$brain_region" {print \$1}' $columns_data};
					if (!$column_num){$columns_num[0] =" "; print "No columns in this category: $time .. $brain_region\n";}
					else {
						chomp $column_num;
						@columns_num = split('\n',$column_num);
						}
					my $median_value;
					print "Looking at $hemis $brain_region\n";
					if ($columns_num[0]eq" "){$median_value = "";}
					else {$median_value = &median(\@tab,\@columns_num);} #&median(@columns_num);}
					print "$median_value\n";
					print $OUT_pd "\t".$median_value;
					}
				}
			}
			print $OUT_pd "\n";
			$gene_count++;
		}
		
	close $OUT_pd;
	my $starttime = localtime();
	print $LOG1 "$starttime: Finished pulling out the exonarray data!!!\n";
	close INFILE;
	qx {rm $temp_outfile};
}

#################################################################################
##  SUBROUTINES for pulling_out_all_data_exonarray ##
#################################################################################

#finding the median value in an array of values
sub median{
	my ($array_nums, $columns) = @_;
	my @array = @$array_nums;
	my @columns_num = @$columns;
#	print "my first value of array is $array[0]\n";
	#print "my first value of array is $columns_num[0]\n";
	my @array_of_values;
#	my $array_num = scalar @array;
	my $array_num = scalar @columns_num;
	print "There are $array_num columns in this group\n";
	print "The last value is $columns_num[$array_num-1]\n";
	my $odd_even;
	my $median;
	my $i;
#	if ($array[0]eq" "){$median = "";print "There is no value\n";}
	if ($columns_num[0]eq" "){$median = "";print "There is no value\n";}
	else{
	for ($i=0;$i<$array_num;$i++){
#		my $current_col = $array[$i];
		my $array_num = $columns_num[$i] - 1;
		my $current_col = "\$"."$array_num";
#		my $current_col = $columns_num[$i];
		#print "my current column is $current_col\n";
#		my $temp = qx {awk '{print $current_col}' $temp_outfile};
		my $temp = $array[$array_num];
		#print $temp."\n";
		#chomp $temp;
		#my @temp_values = split('\n',$temp);
		#if ($gene_list eq "all") {
		#	shift @temp_values;}
		#push (@array_of_values, @temp_values); 
		push (@array_of_values, $temp);
		}
		my $num_values = scalar @array_of_values;
		my @temp_array = sort {$a <=> $b}@array_of_values;
		my $mid = (int($num_values/2));
		print "there are $num_values values to look for the median of\n";
		if ($num_values%2) {
			$odd_even = "odd";
			$median = $temp_array[$mid];
			#print "my median value is $temp_array[$mid] *$median*\n";
			}
		else {
			$odd_even = "even";
			#print "mid is $mid\n";
			$median = (($temp_array[$mid]+$temp_array[$mid-1])/2);
			#print "my median value is $median\n";
			}
		}
	#print "my array has an $odd_even number of values\n";
	return $median;
	}

#================================================================================================
	### DESCRIPTION: To take a file with the median scores of expression data for each brain region/time for a list of genes
	### and provide as options time periods and brain regions of interest and find which genes meet the cutoff for expression and make a list of genes that match
	### Please provide the median score expression data file, the output name for the new file and select the options with values to pull out the genes that meet criteria for which brain regions and time periods of interest
	### EX: <median_scores_expression_file><output_name><STR_R,S1C_L><9,10>  
	### brain regions must include hemispheres ex. STR_R  
	#==================================================================================================\n";	
	
sub criteria_expression_median_genes {
	my ($expression_data, $output_name, $brain_interest, $time_chosen) = @_;
		
	my $outfile_crit_express = "${output_name}_genes.txt";
	print "Making output file $outfile_crit_express.\n";
	$starttime = localtime();
	print $LOG1 "$starttime: Making output file $outfile_crit_express to contain the genes found in the expression data that meet the cutoff of expression in the brain regions and time periods of interest.\n\n";
	open my $OUT_ce, '>', $outfile_crit_express;
	my @genes;

	my @gene_names_list;
	my $final_gene_count;
	my $column;
	my $column_num;
	my @column_nums;
	my @brain_regions;
	my @time_periods;

	my @genders = ('M','F');
	my @hemis = ('L','R');

	@time_periods = (1..15);

	my $cutoff = 6;

	#determine what option to run
	@brain_regions = split(',',$brain_interest);
	my @times = split(',',$time_chosen);
	my @to_look;

	foreach my $time (@times){
		foreach my $region (@brain_regions){
			my $merge = $time."_".$region;
			#print "$merge\n";
			push @to_look, $merge;
			}
		}
	$starttime = localtime();
	print $LOG1 "$starttime: Looking at the following time/brain_regions in the file $expression_data.\n";
	print "Looking at the following time/brain_regions in the file $expression_data.\n";

	my $headers = qx {head -1 $expression_data};
	chomp $headers;
	my @headers = split('\t',$headers); 
	my $total_lines = qx {wc -l < $expression_data};
	chomp $total_lines;

	foreach my $head (@to_look){
		($column) = grep {$headers[$_] =~ /$head/} 0..$#headers; 
		print "The column for $head is $column.\n";
		push @column_nums, $column;
		}

	###Start making file
	print $OUT_ce "Gene\n";

	my $current_line =2; 
	print "Line is $current_line.\n";

	print "To look at each gene in $expression_data and see if any time/region of interest meets the cutoff criteria.\n";
	for (my $i=2;  $i<=$total_lines; $i++){
		my $line = qx {sed -n '${i}p' $expression_data};
		$line =~ s/[\n\r]+//g;
		my @tab = split /\t/,$line;
		my $check = "";
		print "Looking at line $current_line.\n";
		print "Looking at $tab[0] gene currently.\n";
		foreach my $columns (@column_nums){
			if ($check ne "") {next;}
			if ($tab[$columns] >= $cutoff) {
				print $OUT_ce $tab[0]."\n";
				$check = $tab[0];
				}
			}
		$current_line +=1;
		}	
		
	close $OUT_ce;
	$starttime = localtime();
	print $LOG1 "Finished finding the genes that are expressed!!\n\n";
	}
	
	#================================================================================================
	### DESCRIPTION: To make a new file containing only those from a second files list 
	### Please provide the original list of sample and genes, the file with final list of genes to keep, and output file name
	### EX: <genes_list_of_interest> <genes_to_keep> <output_name>
	### OUTPUT: Will provide a new text file containing the sample names and only those corresponding genes from the list to keep
	#==================================================================================================\n";
	
sub extract_gene_list {
	my ($orig_gene_list, $genes_to_keep, $output_name) = @_;
	my $outfile_exonarray = "${output_name}_exonarray_data.txt";
	open my $OUT_E_data, '>', $outfile_exonarray;
	$starttime = localtime();
	print $LOG1 "$starttime: Making file $outfile_exonarray to contain the CNV data for the genes in the file $genes_to_keep.\n";

	my $temp_list = qx {awk '{print \$0}' $genes_to_keep};
	chomp $temp_list;
	my @genes = split('\n',$temp_list);
	$starttime = localtime();
	print $LOG1 "$starttime: This list $genes_to_keep has ",scalar(@genes)," number of genes to keep in it.\n";
	print "This list $genes_to_keep has ",scalar(@genes)," number of genes to keep in it.\n";
	#@genes = grep(s/\s*$//g, @genes);

	my $column;

	open INFILE , '<' , $orig_gene_list or die "Could not open file $orig_gene_list\n$!";

	while (my $line = <INFILE>){
		$line =~ s/[\n\r]+//g;
		if ($line =~ /Gene/) {
			print $OUT_E_data $line."\n"; 
			my @hids = split /\t/,$line; 
			($column) = grep {$hids[$_] =~ /Gene|Genes/} 0..$#hids; 
			#print "The columns containing the gene names is ($column+1).\n";
			next;}
		else {
			my @tab = split /\t/,$line;
			my $numcol = scalar @tab;
			#print "There are $numcol columns.\n";
			my @these_genes = split /,/,$tab[$column];
			my $num_genes = scalar @these_genes;
			#print "This row has $num_genes genes to look at.\n";
			my $new_genes_list="";
			for (my $i=0; $i<$num_genes; $i++){
				$these_genes[$i] =~ s/\"//g;
				#print "Gene $i is **$these_genes[$i]**\n";
				if (grep (/^$these_genes[$i]$/i,@genes)) {
					if ($new_genes_list eq ""){
						$new_genes_list = $these_genes[$i];
						}
					else {
						$new_genes_list=$new_genes_list.",".$these_genes[$i];
						}
					}
				}
			if ($new_genes_list eq ""){next;}
			else {
				my $last_col = ($numcol-1);
				for (my $j=0;$j<$numcol;$j++){
					if ($j == $column && $j == $last_col) {
						#print "adding the genes\n";
						print $OUT_E_data $new_genes_list."\n";
						}
					elsif ($j == $column){
						print $OUT_E_data $new_genes_list."\t";
						}
					elsif ($j != $column && $j == $last_col) {
						print $OUT_E_data $tab[$j]."\n";
						}
					else {
						print $OUT_E_data $tab[$j]."\t";
						}
					}
				}
			}
		}

	close $OUT_E_data;
	close INFILE;
	$starttime = localtime();
	print $LOG1 "$starttime: Finished making the new CNV file.\n\n";
	}