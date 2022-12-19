#!/usr/bin/perl -w
#
# Nov 28, 2000 added option to email results
# 7/4/01 gap option is turned off permanently
#11/11/03 Added SIFT_queue stuff  JGH
# 3/8/04 Don't wait for search to finish, send them format.pl URL  JGH
#-------------------------------------------------------------------------

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.
use DBI;
use Tie::IxHash;
use Digest::MD5 qw (md5 md5_hex md5_base64);
use Getopt::Std;
use File::Basename;
#Added by Jing Hu for timing
use Time::HiRes qw( time );
my $start = time();
#my $SIFT_HOME = "/projects/simnl/SIFT.BACKUP/sift5.0.4" ;# $ENV{'SIFT_HOME'}; # CHANGE
#my $SIFT_HOME = $ENV{'SIFT_HOME'};
my $SIFT_HOME = '/work/NBS/SIFT/SIFTINDEL';
use vars qw($opt_i $opt_c $opt_d $opt_o $opt_p $opt_m $opt_A $opt_B $opt_C $opt_D $opt_E $opt_F $opt_G $opt_H $opt_I $opt_J $opt_K );
getopts("i:c:d:o:p:m:A:B:C:D:E:F:G:H:I:J:K");
my $usage = "usage: 
$0 
        -i <Query indels filename with complete path>
        -c <Human Genome build: hs37 or hs38>
        -d <Human Genome DB directory path>
        -o <output folder (without underscore) - default=<SIFT_HOME>/tmp>
        -p <Optional: prefix name of file>
        -m 1 to output one annotation row per indel: default = 0

        All values should be in local 0 space based coordinates.

        -A 1 to output Gene Name: default: 1
        -B 1 to output Gene Description: default: 1
        -C 1 to output Ensembl Protein Family ID: default:1 
        -D 1 to output Ensembl Protein Family Description: default: 1
        -E 1 to output Ensembl Transcript Status (Known / Novel): default: 1
        -F 1 to output Protein Family Size: default: 1
        -G 1 to output Ka/Ks (Human-mouse): default: 1
        -H 1 to output Ka/Ks (Human-macaque): default: 1
        -I 1 to output OMIM Disease: default: 1
	-J 1 to output 1000 Genomes: default: 1
";
$| = 1;
require "$SIFT_HOME/bin/SIFT_subroutines.pm"; # CHANGE

if(!(
        defined($opt_i) && 
	defined($opt_d) &&
        defined($opt_c))){
    print STDERR $usage;
    die;
}

my $oo1 = defined($opt_A) ? $opt_A:1;
my $oo2 = defined($opt_B) ? $opt_B:1;
my $oo3 = defined($opt_C) ? $opt_C:1;
my $oo4 = defined($opt_D) ? $opt_D:1;
my $oo5 = defined($opt_E) ? $opt_E:1;
my $oo6 = defined($opt_F) ? $opt_F:1;
my $oo7 = defined($opt_G) ? $opt_G:1;
my $oo8 = defined($opt_H) ? $opt_H:1;
my $oo9 = defined($opt_I) ? $opt_I:1;
my $oo10 = defined($opt_J) ? $opt_J:1;
my $output_options = "$oo1,$oo2,$oo3,$oo4,$oo5,$oo6,$oo7,$oo8,$oo9,$oo10";
print "output options $output_options\n";

#       Set file permissions to rw-rw----
#system("umask 006");
my $bin = "$SIFT_HOME/bin"; # CHANGE
#Feb 14, 2013. non-3n indel prediction bin, updated by Jing Hu
my $non3nindelbin	    = "$bin/indel-bin"; # CHANGE
#Feb 4. 2014. 3n indel prediction bin, added by Jing Hu
my $threeNindelbin			="$bin/indel-bin";

my $date = `date`; chomp ($date);
my $key = substr (md5_hex($$ . $date), 0,10);

#my $pid             = $$;
#my $pid = $key;
$pid = "indel" . $$;
if (defined($opt_p)) { $pid = $opt_p; }
my $tmp = defined($opt_o) ? "$opt_o/$pid" : "$SIFT_HOME/tmp/$pid";
#my $tmp = "/mnt3/tmp/$pid"; # TO TAKE ARGUMENTS

#for indel folder with underscore problems (temporary fix)
my $original_tmp = defined($opt_o) ? "$opt_o" : "$SIFT_HOME/tmp";
$tmp =~ s/\_/\-/g;
$parent=$opt_o;
$parent =~ s/\_/\-/g;

print "making directory $tmp\n";
if (! -d $tmp) { system ("mkdir -p $tmp"); }
system ("rm -f $tmp/*"); # remove files already in tmp else it won't work
system ("chmod 777 $tmp");

## SIMNL<2013.01.25>
## Because the input directory file may contain 
## underscores as well, we copy the file to
## the output directory.
my $file_to_use = $tmp . "/" . $$ . ".input";
system("cp $chr_file $file_to_use");
$chr_file = $file_to_use;


my $genome_info_dir = "/mnt1/db"; # TO TAKE ARGUMENTS
my $return_address     = "sift-no-reply\@bii.a-star.edu.sg";
my $snp_classifier     = "perl $bin/Classify_SNPs_Indels_jhu_March2013.pl";
# March 28, 2011, change because it's no longer pid
my $reformat_chrfile   = "perl $bin/reformat_chrfile_032811.pl";
my $detect_indel 	    = "perl $bin/detect_indel.pl";
my $model_transcript   = "perl $bin/model_transcript.pl";
my $small_indels_error_intersect_file = "$tmp/$pid.small_indels_error_intersect.gff";
my $indelfile_to_gff = "perl $bin/indelfile_to_gff.pl";
my $intersect_locations = "java -Xmx1000m -jar $bin/IntersectFeatures.jar ";
my $detect_repeat = "perl $bin/detect_repeat.pl";
#updated by Jing Hu @ Feb 14, 2013
my $non3nindelpred = "python $non3nindelbin/non3nindelPipeline.py";
#Feb 4, 2013 added by Jing Hu
my $threeNindelpred = "python $threeNindelbin/3nindelPipeline.py";
#$foobar = "python $indelbin/foobar.py";
my $genomes1000_intersect_file = "$tmp/$pid.1000Genomes_intersect.gff";
 
# output the beginning text to be used on all pages
#print "Content-type: text/html\n\n";
#print "<body bgcolor=white>\n";

#if ( $ENV{"REQUEST_METHOD"} ne "POST" ) {
#	print "This script should be referenced with a METHOD of POST\n";
#	exit;
#}
=pod
my $QUERY_STRING = "";
read( STDIN, $QUERY_STRING, $ENV{"CONTENT_LENGTH"} );
my %names = &general_parse( $ENV{"CONTENT_TYPE"}, $QUERY_STRING );
my $address;
chomp( $names{address} );
$names{address} =~ s/\s+//;
if ( $names{address} ne "" ) {
	$address = $names{address};
}
=cut
=pod
#get user output options # TO TAKE ARGUMENTS
my $all_transcripts = $names{ALL_TRANSCRIPTS}==1 ? 1 : 0;
my $oo1 = $names{oo1}==1 ? 1 : 0;
my $oo2 = $names{oo2}==1 ? 1 : 0;
my $oo3 = $names{oo3}==1 ? 1 : 0;
my $oo4 = $names{oo4}==1 ? 1 : 0;
my $oo5 = $names{oo5}==1 ? 1 : 0;
my $oo6 = $names{oo6}==1 ? 1 : 0;
my $oo7 = $names{oo7}==1 ? 1 : 0;
my $oo8 = $names{oo8}==1 ? 1 : 0;
my $oo9 = $names{oo9}==1 ? 1 : 0;
my $oo10 = $names{oo10}==1 ? 1:0;
my $output_options = "$oo1,$oo2,$oo3,$oo4,$oo5,$oo6,$oo7,$oo8,$oo9,$oo10";
=cut
## Check for validity of user inputs
## Check for validity of user inputs
my $all_chr_file = $tmp . "/" . $pid . ".chrfile";
my $chr_file = $opt_i; 
#copy($chr_file, $all_chr_file) or die "File $all_chr_file cannot be copied.";
system("cp $chr_file $all_chr_file");
=pod
my $all_chr_file = $tmp . "/$pid.chrfile";
if ( $names{CHR} eq "" && $names{CHR_file} eq "" ) {
	print
"<H1>Error</H1> Please enter some chromosome coordinates with substitutions.<P>\n";
	exit;
}
=cut
=pod
my $organism = $names{organism};
$organism =~ s/^\s+//;
$organism =~ s/\s+$//;
=cut

my $genomebuild = "$opt_c";
$genomebuild =~ s/^\s+//;
my $coding_info_dir = "$opt_d/Coding_info_37_V66";
if(lc($genomebuild) eq "hs38") {$coding_info_dir = "$opt_d/Coding_info_38_V83";}

my $Variation_db_dir = "$opt_d";
my $bin_file = "$coding_info_dir/bins.list"; 
my $transcript_db = "$Variation_db_dir/TRANSCRIPT_DB.sqlite";
=pod
if ($organism eq "hs36") {
    $bin_file = $genome_info_dir . "/Coding_info_36/bins.list"; 
    $transcript_db = $genome_info_dir . "Human_db_37/TRANSCRIPT_DB.sqlite";
} elsif ($organism eq "hs37") {
	$bin_file = $genome_info_dir . "/Coding_info_37_V66/bins.list";
#	$transcript_db = $genome_info_dir . "/Human_db_37/TRANSCRIPT_DB.sqlite";
   #gina changed for V66 DB 
	$transcript_db = $genome_info_dir . "/TRANSCRIPT_DB_V66.sqlite";
}
=cut
=pod
if ( $organism =~ /Select Organism/i ) {
	print
"<H1>Error</H1> Please select organism after pressing back button on your browser.<P>\n";
	exit;
}
=cut
=pod
open( CHR_FILE, ">$all_chr_file" );
my $all_chr_string;
if ( $names{CHR_file} ne "" ) {
	$names{CHR_file} =~ s/\r/\n/g;
	$names{CHR_file} =~ tr/A-Z/a-z/;
	if ( $names{CHR_file} =~ /\d/ ) {
		print CHR_FILE uc( $names{CHR_file} );
	}
	$all_chr_string = uc( $names{CHR_file} ), "\t";
}
if ( $names{CHR} ne "" ) {
	$names{CHR} =~ tr/A-Z/a-z/;
	$names{CHR} =~ s/^\n//;
	$names{CHR} =~ s/\r/\n/g;
	print CHR_FILE uc( $names{CHR} );
	$all_chr_string .= uc( $names{CHR} ), "\t";
}
close(CHR_FILE);
chmod 0777, $all_chr_file;
=cut
### 
### If user input is > 1000, the website throws out a message that the file is too large and to do 
### (1) Intersect With Coding webpage before submitting to the SIFT_chr_coords_submit. 
###
my $number_of_inputs = 0;
open(CHR_FILE, "<$all_chr_file");
while (my $line = <CHR_FILE>) {
    chomp $line;
    $line =~ s/^\s+//;
    $line =~ s/\s+$//;
    if ($line ne "" && $line !~ /ONTENT/) {
	$number_of_inputs++;
    }
}
close(CHR_FILE);

#####gina ADD preprocessor########
system("$bin/indels-prep.py $all_chr_file");
my $prepRet=$?;
system("mv $all_chr_file $all_chr_file".".org");
system("cp $all_chr_file".".ok"." $all_chr_file");
#####gina ADD preprocessor END########
=pod
if ($number_of_inputs > 1000) {
    print "The file you provided has more than 1000 entries.<BR>";
    print "Please partition into separate files of not more than 1000 each and submit separately.<BR>";
    print "<a href=\"/www/SIFT_chr_coords_indels_submit.html\">Back</a><P>";
    exit(0);
}
=cut

#check input validity and reformat chrfile
chmod 0777, $all_chr_file;
#check input validity and reformat chrfile
system("$reformat_chrfile $all_chr_file $tmp");
print "$reformat_chrfile $all_chr_file $tmp\n";
#system("$reformat_chrfile $all_chr_file");

select(STDOUT);
$|++;

#display wait message
$| = 1;
#print
#"Your job id is $pid and is currently running.\nIf your browser times out before results are shown, your results will be stored at $url_address/www/sift/tmp/$pid/$pid\_predictions.html </font>" . "  for 24 hours." . 
#"<BR><BR>
#Problems? Please contact <A HREF=\"contact.pl\">us<A> with your job id.\n";
#select(STDOUT);


#convert chrfile to gff for intersecting with small indel errors file.
system("$indelfile_to_gff $all_chr_file");
my $all_chr_gff_file = "$all_chr_file.gff";

#Now intersect gff chrfile with small indels ref file


my $small_indels_error_gff = "$opt_d/small_indel_errors_ref.gff";
my $genomes1000_gff = "$opt_d/all_indels_bed_10bp.gff";
=pod
if ($organism eq "hs36") {
    $small_indels_error_gff =  $genome_info_dir . "/Human_db_36/small_indel_errors_ref.gff";
	$Human_db_dir = $genome_info_dir . "/Human_db_36/"; 
	$genomes1000_gff = $genome_info_dir . "/Human_db_36/all_indels_bed36_10bp.gff";

} elsif ($organism eq "hs37") {
	$small_indels_error_gff = $genome_info_dir . "/Human_db_37/small_indel_errors_ref.gff";
	$Human_db_dir = $genome_info_dir . "/Human_db_37/";
	$genomes1000_gff  = $genome_info_dir . "/Human_db_37/all_indels_bed37_10bp.gff";

}
=cut
system("$intersect_locations $all_chr_gff_file gff $small_indels_error_gff gff simple > $small_indels_error_intersect_file");
system ("$intersect_locations $all_chr_gff_file gff $genomes1000_gff gff simple > $genomes1000_intersect_file");

#create an index %user_small_indels_error_index with key as "chr#startindex"  stop index redundant also problem with Intersect script not detecting intersect when start = stop i.e. insertion.
my %user_small_indels_error_index = ();
if (-e $small_indels_error_intersect_file){
	open (SMALL_INDELS,"$small_indels_error_intersect_file"); 
	while (<SMALL_INDELS>){
		next if ($_ =~ /^\d+|^Total/i);
		chomp;
		my @elts = split /\t/, $_;
		my $chr = $elts[1];
		my $start = $elts[2];
		my $key = "$chr\t$start";
		$user_small_indels_error_index{"$key"} = 1;	
	}
	close (SMALL_INDELS);
}
my %user_1000G_index = ();
if (-e $genomes1000_intersect_file) {
	open (GENOMES1000, "$genomes1000_intersect_file");
	while (<GENOMES1000>) {
		next if ($_ =~ /^\d+|^Total/i);
                chomp;
                my @elts = split /\t/, $_;
                my $chr = $elts[1];                my $start = $elts[2];
                my $key = "$chr\t$start";
                $user_1000G_index{"$key"} = $elts[4] . "." . $elts[11];
	}
	close (GENOMES1000);
}

#Creat binned SNP files for SNP Classifier
print "MAP COORDS: perl $bin/map_coords_to_bin_indels.pl $bin_file $all_chr_file $tmp $pid\n";
system("perl $bin/map_coords_to_bin_indels.pl $bin_file $all_chr_file $tmp $pid");

# At this point, snp_chr_map_file.txt is created, snps_chrNum_binNum_binStart_binEnd.txt files are created
# Each of the text files may contain one or more line:
# 1       63760421        63760421        1       GCGGCT  #USER COMMENT 1
# 2       11205130        11205130        1       ACACACACACAC

#create index of old coords to new coords using pid_old_new_coords_map.txt
print "TRY THIS: $tmp/$pid\_old_new_coords_map.txt\n";
open( COORDS_MAP, "$tmp/$pid\_old_new_coords_map.txt" )
  || die("Cannot open old-new coords map file");
while (<COORDS_MAP>) {
	chomp;
	@elts       = split /\t/, $_;
	$old_coords = $elts[0];
	$new_coords = $elts[1];
	@elts2 = split /,/, $new_coords;
	$chr = $elts2[0];
	$new_start = $elts2[1];
	$new_stop = $elts2[2];
	###changed by Jing Hu to fix the bug 
	# if there are two indels : 1,171174749,171174749,-1,AGTTCGA, and 1,171174749,171174749,1,AAT,
	# previous program cannot tell the difference.
	$newallele = $elts2[4];
	$key = "$chr\t$new_start\t$new_stop\t$newallele";
	$index_old_new_coords{$key} = $old_coords;
} #end while
close(COORDS_MAP);

#Run snp classifier on files in snp_chr_map_file
open( SNP_CHR_MAP_FILE, "$tmp/$pid\_snp_chr_map_file.txt" )
  || die("Cannot open SNP_CHR_MAP_FILE");
my $full_coding_info_dir = $coding_info_dir;
#if ($organism eq "hs36") {
#	$full_coding_info_dir = $genome_info_dir . "/Coding_info_36";
#} elsif ($organism eq "hs37") {
#	$full_coding_info_dir = $genome_info_dir . "/Coding_info_37_V66";
#}

# SIM to prevent timeout
#print "<BR>";
#print "<b>Please do not refresh or close this window.</b>  Periods will be occassionally printed below to prevent browser time-out.<BR><BR>";
#print "<BR>Working..";

while (<SNP_CHR_MAP_FILE>) {
	chomp;
	@cols             = split /\t/, $_;
	$chromosome       = $cols[0];
	$bucket           = $cols[1];
	$snp_file         = "$tmp/$cols[2]";
	$coding_info_file = "$full_coding_info_dir/$cols[3]";
	$chr_fasta_file   = "$full_coding_info_dir/$cols[4]";

	# SIM TO PREVENT TIMEOUT
#	print ".";

# Pauline this commented out Sept 9 2010 but need to put it back in to make it faster
#	 $coding_info_file = "$coding_info_dir/$organism/chr$chromosome/$cols[3]";        $chr_fasta_file   = "$coding_info_dir/$organism/chr$chromosome/$cols[4]";        #print "$snp_classifier -s $snp_file -c $coding_info_file -n $chr_fasta_file -o $tmp/$pid\_snps_classified_chr$chromosome\_bin$bucket<br>";


	#print "$snp_classifier -s $snp_file -c $coding_info_file -n $chr_fasta_file -o $tmp/$pid\_snps_classified_chr$chromosome\_bin$bucket<br>";
print "$snp_classifier -s $snp_file -c $coding_info_file -n $chr_fasta_file -o $tmp/$pid\_snps_classified_chr$chromosome\_bin$bucket\n"; 
system ("$snp_classifier -s $snp_file -c $coding_info_file -n $chr_fasta_file -o $tmp/$pid\_snps_classified_chr$chromosome\_bin$bucket"); 


}
close(SNP_CHR_MAP_FILE);




#Parse denormal files to extract coding info for input coords
open( SNP_CHR_MAP_FILE, "$tmp/$pid\_snp_chr_map_file.txt" )
  || die("Cannot open SNP_CHR_MAP_FILE");
#Modified by Jing Hu @ Feb 4, 2013. ####################
my $chr_file_flanking = $tmp . "/$pid.chrfile.flanking";
open( CHR_FILE_FLANKING, ">$chr_file_flanking" );
#########################################################
my %coord_result_hash;
my $en_trancript_id_string;
while (<SNP_CHR_MAP_FILE>) {
	chomp;
	@cols             = split /\t/, $_;
	$chromosome       = $cols[0];
	$bucket           = $cols[1];
	$snp_file         = "$tmp/$cols[2]";
# Pauline modified Sept 2010 because it's not partitioned 
#	$coding_info_file = "$coding_info_dir/$organism/chr$chromosome/$cols[3]";
#	$chr_fasta_file   = "$coding_info_dir/$organism/chr$chromosome/$cols[4]";
#	$coding_info_file = "$coding_info_dir/$organism/$cols[3]";
#        $chr_fasta_file   = "$coding_info_dir/$organism/$cols[4]";


	$denormal_file    =
	  "$tmp/$pid\_snps_classified_chr$chromosome\_bin$bucket.denormal";
#	$normal_file    =
#          "$tmp/$pid\_snps_classified_chr$chromosome\_bin$bucket.normal";
#### 2013.07.03: For some reason, $normal_file was declared but never used. This throws out warnings when running the script          

	open( DENORMAL, $denormal_file ) || die "(denormal file does not exist)";
	while (<DENORMAL>) {
		@elts        = split /\t/, $_;
		$coord       = $elts[1];
		@coord_elts  = split /\:/, $coord;
		$allele      = $coord_elts[scalar @coord_elts - 1];
		$coord       = join( ',', @coord_elts );
		@coord_elts2 = split( '-', $coord, 2 );
		$coord       = join( ',', @coord_elts2 );
		@elts2 = split /,/, $coord;
		$new_start = $elts2[0];
		$new_stop = $elts2[1];
		#changed by Jing Hu to fix the previous problem, which cannot handle when two indel has the same coordinates
		$direction = $elts2[2];
		#$allele = $elts2[-1];
		$allele = $elts2[3];
		$new_allele = "";
		if ($direction eq "-1" && $new_start==$new_stop){
			my $i;
			for ($i = length($allele)-1 ; $i >= 0; $i--){ 
				if (substr($allele, $i, 1) eq "A") {
				$new_allele = $new_allele."T";
				}
			elsif (substr($allele, $i, 1) eq "T") {
				$new_allele = $new_allele."A";
				}
			elsif (substr($allele, $i, 1) eq "C") {
				$new_allele = $new_allele."G";
				}
			elsif (substr($allele, $i, 1) eq "G") {
				$new_allele = $new_allele."C";
				}
			}
		}
		elsif ($new_start!=$new_stop && ($allele eq "-" || $allele eq "-/")){
			$new_allele = "-/";
			}
		else{
			$new_allele = $allele;
			}
				
		$key = "$chromosome\t$new_start\t$new_stop\t$new_allele";
		$old_coord= $index_old_new_coords{$key};
		$coord       = $old_coord;
		$nu_change = $elts[2];
		

                if ($nu_change =~ /([ATGC]+) ([ATGC]+) ([ATGC]+)/){
                        $left_flank = uc substr ($1,-5,5);
                        $deletion = lc $2;
                        $right_flank = uc substr($3,0,5);
                        $nu_change_modified = "$left_flank-$deletion-$right_flank";
                        
                        #Added by Jing Hu @ 2/4/2013 to generate flanking sequences of 50 bases on each side
                        $long_left_flank = uc substr ($1,-50,50);
                        $deletion = lc $2;
                        $long_right_flank = uc substr($3,0,50);
                        $long_nu_change_modified = "$long_left_flank-$deletion-$long_right_flank";
                }
                elsif($nu_change =~ /([ATGC]+) - ([ATGC]+)/){
                        $left_flank = lc substr ($1,-5,5);
                        $insertion = uc $allele;
                        $right_flank = lc substr($2,0,5);
                        $nu_change_modified = "$left_flank-$insertion-$right_flank";
                        
                        #Added by Jing Hu @ 2/4/2013 to generate flanking sequences of 50 bases on each side
                        $long_left_flank = lc substr ($1,-50,50);
                        $insertion = uc $allele;
                        $long_right_flank = lc substr($2,0,50);
                        $long_nu_change_modified = "$long_left_flank-$insertion-$long_right_flank";
                        
                }
		$nu_change = $nu_change_modified;
		#added by Jing HU @ 2/4/2013
		print CHR_FILE_FLANKING "$coord\t$long_nu_change_modified\n";
		
		# SIM TEST
		#print "After modification: $nu_change<BR>";
		# END SIM TEST

		if ( $_ =~ /AA_DELETION | AA_INSERTION | EARLYSTOP | SYNONYMOUS | FRAMESHIFT/x ) {  #Updated by Jing Hu @ March 2013.
			$subst = $&;
			$region = "CDS";
                        if ( $_ =~ /^(\d+).+?ENSG(\d+).+?ENST(\d+).+/) {
				$hit_id = $1;
                                $en_transcript_id = "ENST$3";
				$ensg = "ENSG$2";
                                $result = "$hit_id\t$ensg\t$en_transcript_id\t$subst\t$region";
				#append seq (original and modified) to result
				$original_seq = "$tmp/$pid\_$en_transcript_id.fa";
				$modified_seq = "$tmp/$pid\_$en_transcript_id\_$hit_id.fa";
				if (-e $modified_seq && -e $original_seq){
					$result.="\t$original_seq\t$modified_seq\t$nu_change";
				}
				push @{ $coord_result_hash{$coord} }, $result;		#hash of arrays: multiple ENSTs per coord
                        }
                }
		elsif ( $_ =~ /INTRON|PROMOTOR|DOWNSTREAM|UPSTREAM/ix ) {
                        $subst = "";
			$region = "$&";
                        if ( $_ =~ /^(\d+).+?ENSG(\d+).+?ENST(\d+).+/) {
                                $hit_id = $1;
                                $en_transcript_id = "ENST$3";
                                $ensg = "ENSG$2";
                                $result = "$hit_id\t$ensg\t$en_transcript_id\t$subst\t$region\t\t\t$nu_change";
                                #append seq (original and modified) to result
                                push @{ $coord_result_hash{$coord} }, $result;          #hash of arrays: multiple ENSTs per coord
                        }
                }
		
		elsif ( $_ =~ /CDS/i ) {
                        $subst = "";
                        $region = "CDS";
                        if ( $_ =~ /^(\d+).+?ENSG(\d+).+?ENST(\d+).+/) {
                                $hit_id = $1;
                                $en_transcript_id = "ENST$3";
                                $ensg = "ENSG$2";
                                $result = "$hit_id\t$ensg\t$en_transcript_id\t$subst\t$region\t\t\t$nu_change";
                                #append seq (original and modified) to result
                                push @{ $coord_result_hash{$coord} }, $result;          #hash of arrays: multiple ENSTs per coord
                        }
                }
        
			
	}
	close (DENORMAL);


	#unlink ("$denormal_file");
	#unlink ("$normal_file");
}

close CHR_FILE_FLANKING;

# SIM TEST
#print "<BR>";
#foreach my $key (keys %coord_result_hash) {
#    my $arrRef = $coord_result_hash{$key};
#    my @array = @{$arrRef};
#    print "KEY: $key<BR>";
#    foreach my $result (@array) {
#	print "   $result<BR>";
#    }
#    print "<P>";
#}
# END TEST SIM


#Modified by Jing Hu @ Feb 21, 2013. ####################
#This file keeps these 3n indels which introduce stop codons
my $chr_file_3nSTOP = $tmp . "/$pid.chrfile.3nSTOP";
open(CHR_FILE_3nSTOP, ">$chr_file_3nSTOP" );

#create original - modified protein files
my $repeat_detect_input_file = "$tmp/$pid\_repeat_detect_input_file.txt";
open (REPEAT_INPUT,">$repeat_detect_input_file");

foreach $coord (keys %coord_result_hash){
	my @elts;
	my $repeat_input_row = "";
	my @indel_result;
	my @indel_location_features;
	my $chr = (split /,/, $coord)[0];
	my $indel_start = (split /,/, $coord)[1];
	my $indel_stop = (split /,/, $coord)[2];
	my $allele = (split /,/, $coord)[4];
	@result_set = @{ $coord_result_hash{$coord} };
	foreach $result (@result_set){
		chomp $result;
		@elts = split /\t/,$result; #print "RESULT: $result\n";
		my $hit_id = $elts[0];
		my $ensg = $elts[1];
		my $enst = $elts[2];
		my $subst = $elts[3];
		my $region = $elts[4];
		my $original_seq = $elts[5];
		my $modified_seq = $elts[6];
#		print "MODEL_TRANSCRIPT: $model_transcript $enst $indel_start $indel_stop $transcript_db\n";
#		@indel_location_features = split /\t/, `$model_transcript $enst $indel_start $indel_stop $transcript_db`;
		@indel_location_features = split /\t/, `$model_transcript $enst $indel_start $indel_stop $transcript_db`;


		### TESTING
#		foreach my $test (@indel_location_features) {
#		    chomp $test;
#		    print "TEST: $test\n";
#		}
		### END TESTING



		$indel_start_CDS = $indel_location_features[0];
		$indel_stop_CDS = $indel_location_features[1];
		$indel_loc_pct =  $indel_location_features[2];
		$indel_span =  $indel_location_features[3];
#		$num_5UTR = $indel_location_features[4];
		$num_exon = $indel_location_features[5];
#		$num_intron = $indel_location_features[6];
#		$num_3UTR = $indel_location_features[7];
#		$tx_length = $indel_location_features[8];
#		$CDS_length = $indel_location_features[9];
#		$protein_length = $indel_location_features[10];
		$last_exon_length = $indel_location_features[11];
		$visual_tx = $indel_location_features[12];
		$NMD_protein_region = int(($last_exon_length + 50)/3);
		$result.="\t$indel_start_CDS\t$indel_stop_CDS\t$indel_loc_pct\t$indel_span\t$visual_tx";
		my $outfile = "";
		if (-e $original_seq && -e $modified_seq){
		    $outfile = "$tmp/$pid\_$enst\_$hit_id\_comparison.fasta"; #print "$original_seq $modified_seq\n";
			@indel_result = split /\t/, `$detect_indel $original_seq $modified_seq $outfile`;
			$aa_coords_original =  $indel_result[0];
			$aa_coords_modified =  $indel_result[1];
			$aa_change_original =  $indel_result[2];
			$aa_change_modified =  $indel_result[3];
			$original_length = $indel_result[4];
			$modified_length = $indel_result[5];
			#updated by Jing Hu @ Feb 21, 2013 to output 3n indels which introduce early stop codon
			my $diffInd = $indel_stop - $indel_start;
			my $lenAllele = length $allele;
			#if it is 3n indel, we determine if the size change in protein sequence is higher
			#than simple AA insertion/deletion
			if ($diffInd!=0 && $diffInd % 3 == 0 || $diffInd==0 && $lenAllele % 3 == 0) {
				my $sizeChange = $diffInd;
				if ($diffInd==0){
					$sizeChange = $lenAllele;
					}
				$sizeChange = $sizeChange / 3; # this is the possible AA size change due to simple AA insertion/deletion
				if ($modified_length) { chomp $modified_length; }
#				print "ORIG VS MOD: $original_length, $modified_length\n";
				my $AAsizeChange =  $original_length -$modified_length;  #this is the actual AA size change
				
				if ($AAsizeChange > 0 and $AAsizeChange>$sizeChange){
					print CHR_FILE_3nSTOP "$coord\tEARLYSTOP\n";
				}
			}
			
			if ($num_exon < 2) {
				$causes_NMD = "N\/A";
			}
			elsif ($original_length - $modified_length > $NMD_protein_region){
				$causes_NMD = "YES";
			}
			else{
				$causes_NMD = "NO";
			}
#			chomp $aa_coords_original,$aa_coords_modified,$aa_change_original,$aa_change_modified;	
			chomp $aa_coords_original; chomp $aa_coords_modified; chomp $aa_change_original; chomp $aa_change_modified;	
#			chomp $aa_coords_original,$aa_coords_modified,$aa_change_original,$aa_change_modified;	
			$result.="\t$aa_coords_original\t$aa_change_original->$aa_change_modified\t$causes_NMD";
			$repeat_input_row = "$chr\t$enst\t$indel_start\t$indel_stop\t".join ("\t", split(/\-/,$aa_coords_original) )."\t$original_seq\t$modified_seq";
			print REPEAT_INPUT "$repeat_input_row\n";	
		}
	}	
	$coord_result_hash{$coord} = [ @result_set ];
	
}
close (REPEAT_INPUT);
close (CHR_FILE_3nSTOP);

#Run the detect_repeat script on the repeat_input file created in previous loop
my $repeat_detect_output_file = "$tmp/$pid\_repeat_detect_output_file.txt";
#system("$detect_repeat $repeat_detect_input_file $Human_db_dir/repeat_db.gff > $repeat_detect_output_file");
system("$detect_repeat $repeat_detect_input_file $Variation_db_dir/repeat_db.gff > $repeat_detect_output_file");

my %repeat_index;
if (-e $repeat_detect_output_file){
	open (REPEAT_OUTPUT,"$repeat_detect_output_file");
	while (<REPEAT_OUTPUT>){
		chomp;
		my @elts = split /\t/, $_;
		my $chr = $elts[0];
		my $indel_start  = $elts[2];
		if ($elts[10] ne "" && $elts[11] ne ""){	
			my $repeat_change = $elts[10]."-->".$elts[11];
			$repeat_index{"$chr\t$indel_start"} = $repeat_change;
		}
	}	
}


## Check that this IP address hasn't been used too much
#my $IP_address = $ENV{REMOTE_ADDR};

#       print "<HR>" . $IP_address . "<BR></HR> ";
#my $remote_host = $ENV{REMOTE_HOST};
=pod
my $ip_counts =
`cat  /home/blocks/apache/logs/access_log  | grep POST | grep $IP_address | wc -l `;
chomp($ip_counts);
if ( $ip_counts == "" ) {
	$ip_counts =
`cat /home/blocks/apache/logs/access_log  | grep POST | grep $remote_host | wc -l`;
	chomp($ip_counts);
}
=cut
#       print $ip_counts. "<BR>";
#my $upper_limit = 50;
#if ($address && $address ne "" && $address =~ /.*\@.*/) {
#	$upper_limit = 1000;
#}
#if ( $ip_counts > $upper_limit ) {
#	print "<H1>Your computer has exceeded its daily limit.</H1><BR>";
#	print
#"Please download <A HREF=\"/\">SIFT software</A HREF> directly to your computer or <A HREF=\"/sift-bin/contact.pl\">contact</A HREF> us so that we can help you.  Thank you for using SIFT. <BR>";
#	exit;
#}

## SIFT prediction operations

#$exp_option = 1;

#$info = $names{info};
#$comments            = "$tmp/$pid.comments";
# Pauline commented out for increased security Feb 17 2011 
#$out                 = $tmp . "/$$.siftresults";
#$out = $tmp . "/$key.siftresults";



#####################################################################################################
# This code block was added by Jing Hu to parse the indel prediction results and embed results
# into a hash table
# Added by Jing Hu @ Franklin & Marshall College, June 9, 2011.
#####################################################################################################
#added by Jing HU
#system("$non3nindelpred $all_chr_file.orginal $tmp &");
my $subfolder = time(); $subfolder =~ s/(.*)\..*/$1/;
print "subfolder: $subfolder\n";
my $tmp_subfolder= "$tmp/$subfolder";
print "mkdir -p $tmp_subfolder\n";
system("mkdir -p $tmp_subfolder");
system("$non3nindelpred $all_chr_file.orginal $tmp $non3nindelbin $coding_info_dir $tmp_subfolder $genomebuild");
#Feb 4, 2013, added by Jing Hu for 3nindel prediction
#system("$threeNindelpred $chr_file_flanking  $tmp  &");
my $threeNresults = system("$threeNindelpred $chr_file_flanking $tmp $threeNindelbin $coding_info_dir $genomebuild");
#print "THREE N RESULTS: $threeNresults\n";

my %indelClassificationResults;
#read prediction results files
#my $start = time();  #start the timer
my $periods_printed = 0;


=pod
while (not (-o "$all_chr_file.orginal.predictions" or -o "$all_chr_file.orginal.nonExonIndels" or -o "$all_chr_file.orginal.missingIncompleteTranscripts" ) ) {
	
	print(".");
	$periods_printed++;
	sleep(2);
	
	}

#added on Feb 4, 2013 by Jing Hu to wait for 3n indel prediction results
while (not -o "$all_chr_file.3n.predictions" ) {
	
	print(".");
	$periods_printed++;
	sleep(2);
	}

print "<BR><BR>";

print "$all_chr_file.3n.predictions\n";
print "$all_chr_file.orginal.predictions\n";
print "$all_chr_file.orginal.nonExonIndels\n";
print "$all_chr_file.orginal.missingIncompleteTranscripts\n";
=cut

#added by Jing Hu @ Feb 4, 2013 to read results of predictions for 3n indels
if (-o "$all_chr_file.3n.predictions" ){
	open( PREDFILE, "$all_chr_file.3n.predictions" );
	
	while (<PREDFILE>){
		chomp;
		my @predictions = split /\s+/, $_;
		my $indel = $predictions[0];
		
		#3n indel prediction results have indels as chrNum,startIdx,endIdx,orientation,allele
		my @indelParts = split /,/, $indel;
		my $newindelformat;
		if ( $indelParts[4] eq "/"){
			$newindelformat = $indelParts[0].",".$indelParts[1].",".$indelParts[2].",1,-".$indelParts[4];
			}
		else{
			$newindelformat = $indelParts[0].",".$indelParts[1].",".$indelParts[2].",1,".$indelParts[4];
			}
		my $affectedEnsgID = $predictions[1];
		my $effect  = $predictions[2];
		my $confidence = $predictions[3];
		my $classificationRule = $predictions[4];
		my @array = ($effect, $confidence, $classificationRule); 
		@{$indelClassificationResults{$newindelformat}} = @array;
	}
	close(PREDFILE);
	}

	
if (-o "$all_chr_file.orginal.predictions" ){
	open( PREDFILE, "$all_chr_file.orginal.predictions" );
	
	while (<PREDFILE>){
		chomp;
		my @predictions = split /\s+/, $_;
		my $indel = $predictions[0];
		
		#non-3n indel prediction results have indels as chrNum,startIdx,endIdx,allele
		my @indelParts = split /,/, $indel;
		my $newindelformat;
		if ( $indelParts[3] eq "/"){
			$newindelformat = $indelParts[0].",".$indelParts[1].",".$indelParts[2].",1,-".$indelParts[3];
			}
		else{
			$newindelformat = $indelParts[0].",".$indelParts[1].",".$indelParts[2].",1,".$indelParts[3];
			}
		my $affectedEnsgID = $predictions[1];
		my $effect  = $predictions[2];
		my $confidence = $predictions[3];
		my $classificationRule = $predictions[4];
		my @array = ($effect, $confidence, $classificationRule); 
		@{$indelClassificationResults{$newindelformat}} = @array;
	}
	close(PREDFILE);
	}
	

if (-o "$all_chr_file.orginal.nonExonIndels" ){
	open( NONEXONFILE, "$all_chr_file.orginal.nonExonIndels" );
	while (<NONEXONFILE>){
		chomp;
		my $indel = $_;
		
		my @indelParts = split /,/, $indel;
		my $newindelformat;
		if ( $indelParts[3] eq "/"){
			$newindelformat = $indelParts[0].",".$indelParts[1].",".$indelParts[2].",1,-".$indelParts[3];
			}
		else{
			$newindelformat = $indelParts[0].",".$indelParts[1].",".$indelParts[2].",1,".$indelParts[3];
			}
		my @array = ("NONEXON", "NA", "NA"); 
		@{$indelClassificationResults{$newindelformat}} = @array;
	}
	close(NONEXONFILE);
}

if (-o "$all_chr_file.orginal.missingIncompleteTranscripts" ){
	open( MISSFILE, "$all_chr_file.orginal.missingIncompleteTranscripts" );
	while (<MISSFILE>){
		chomp;
		my $indel = $_;
		
		my @indelParts = split /,/, $indel;
		my $newindelformat;
		if ( $indelParts[3] eq "/"){
			$newindelformat = $indelParts[0].",".$indelParts[1].",".$indelParts[2].",1,-".$indelParts[3];
			}
		else{
			$newindelformat = $indelParts[0].",".$indelParts[1].",".$indelParts[2].",1,".$indelParts[3];
			}
	
		my @array = ("MISSING_INCOMPLETE_TRANSCRIOPTS", "NA", "NA"); 
		@{$indelClassificationResults{$newindelformat}} = @array;
	}
	close(MISSFILE);
}



######  Calling the program #########

#print "<A NAME=top><H1><center>S<font color=red>I</font>F<font
#color=blue>T</font> results</center></H1></A><BR>\n";
#if ( $address ne "" ) {
#	print "Results will also be mailed to $address.<BR><BR>\n";
#}

select(STDOUT);
$|       = 1;
#$counter = -1;

#print "<BR><BR>";

#build combined classification (no predictions yet for indels)  file from  %coord_result_hash, hash of arrays. Also get gene infor from Human_supp db
#$db_supp = DBI->connect( "dbi:SQLite:dbname=$Human_db_dir/Human_Supp.sqlite","", "", { RaiseError => 1, AutoCommit => 1 } );

$db_supp = DBI->connect( "dbi:SQLite:dbname=$Variation_db_dir/Human_Supp.sqlite","", "", { RaiseError => 1, AutoCommit => 1 } );
$db_supp->do('PRAGMA synchronous=1');
$sth_db_supp_geneinfo = $db_supp->prepare("select * from GENE_INFO where ENST = ?");

$combined_siftresults_file = "$pid.combined_siftresults"; print "$combined_siftresults_file\n";
open (CLASSIFICATION_FILE,">$tmp/$pid.classification_file");
foreach $coord (keys %coord_result_hash) {
	#added by Jing Hu, June 13, 2011.#########################################
	my @elements = split /,/, $coord;
	my $chrNum = $elements[0];
	my $startIdx = $elements[1];
	my $endIdx = $elements[2];
	my $direction = $elements[3];
	my $allele = $elements[4];

	my $newallele = "";
	if ($direction eq "-1" && $startIdx==$endIdx){
		my $i;
		#Fix the bug of reverse: Jing Hu @ 2/3/2013
		#for ($i = 0 ; $i < length($allele); $i++){ 
		for ($i = length($allele)-1 ; $i >= 0; $i--){
			if (substr($allele, $i, 1) eq "A") {
				$newallele = $newallele."T";
				}
			elsif (substr($allele, $i, 1) eq "T") {
				$newallele = $newallele."A";
				}
			elsif (substr($allele, $i, 1) eq "C") {
				$newallele = $newallele."G";
				}
			elsif (substr($allele, $i, 1) eq "G") {
				$newallele = $newallele."C";
				}
			}
		}
	else{
		$newallele = $allele;
		}

	my $newindel; #indel format in positive chain
	$newindel = $chrNum.",".$startIdx.",".$endIdx.","."1".",".$newallele;
		
	my $pred = $indelClassificationResults{$newindel}[0];
	if (!$pred) { $pred = ""; }
	if ($pred eq "NONEXON" or $pred eq "MISSING_INCOMPLETE_TRANSCRIOPTS"){
		$pred = "NA";
		}
	elsif ($pred eq "disease"){
		$pred = "damaging";
		}
	my $confidence = $indelClassificationResults{$newindel}[1];
	if (!$confidence) { $confidence = ""; }
	my $classificationPath = $indelClassificationResults{$newindel}[2];
	if ($classificationPath && $classificationPath ne "NA" &&  $classificationPath ne ""){ # Need to check if initialized (eg. introns do not have ClassificationPath)
		my @alist = split /_/, $classificationPath;
		$classificationPath = $alist[0]."_CP".$alist[3];
	} if (!$classificationPath) { $classificationPath = ""; } 
	##########################################################################

	if ($subst eq "EARLYSTOP") {
	    print "EARLY STOP: $pred, $confidence, $classificationPath\n";
	}


	my @coord_elts = split /,/, $coord;
	
	#warning if user coordinates intersect with small indel errors
	my $chr = $coord_elts[0];
	my $coord_start = $coord_elts[1];
	my $key_for_small_indels = "$chr\t$coord_start";
	if ($user_small_indels_error_index{$key_for_small_indels} && $user_small_indels_error_index{$key_for_small_indels} == 1){
		$small_indel_warning = "Warning: NCBI reference miscall!";
	}
	else{
		$small_indel_warning = "";
	}

	if (exists ($user_1000G_index{$key_for_small_indels})) {
		$genome1000_annot = $user_1000G_index{$key_for_small_indels};
	} else {
		$genome1000_annot = "";
	}
	
	my $repeat_change = $repeat_index{$key_for_small_indels}; if (!$repeat_change) { $repeat_change = ""; }
	#Populate CLASSIFICATION_FILE with results obtained (multiple results per coord for multiple ENSTs hit)
	@result_set = @{ $coord_result_hash{$coord} };
	foreach $result (@result_set){
		chomp $result;
		@elts             = split /\t/, $result;
		$comparison_file = "";
		$hit_id = $elts[0]; if (!$hit_id) { $hit_id = ""; }
		$ensg = $elts[1];if (!$ensg) { $ensg = ""; }
		$enst = $elts[2];if (!$enst) { $enst = ""; }
		$subst = $elts[3];if (!$subst) { $subst = ""; }
		$region = $elts[4]; if (!$region) { $region = ""; }
		$original_seq = $elts[5]; if (!$original_seq) { $original_seq = ""; }
		$modified_seq = $elts[6]; if (!$modified_seq) { $modified_seq = ""; }
		$nu_change = $elts[7]; if (!$nu_change) { $nu_change = ""; }
		$indel_start_CDS = $elts[8];if (!$indel_start_CDS) { $indel_start_CDS = ""; }
		$indel_stop_CDS = $elts[9]; if (!$indel_stop_CDS) { $indel_stop_CDS = ""; }
		$indel_loc_pct = $elts[10]; if (!$indel_loc_pct) { $indel_loc_pct = ""; }
		$indel_span = $elts[11]; if (!$indel_span) { $indel_span = ""; }
		$visual_tx = $elts[12]; if (!$visual_tx) { $visual_tx = ""; }
		$coords_change = $elts[13]; if (!$coords_change) { $coords_change = ""; }
		$aa_change = $elts[14]; if (!$aa_change) { $aa_change = ""; }
		$causes_NMD = $elts[15];  if (!$causes_NMD) { $causes_NMD = ""; }
		$sth_db_supp_geneinfo->execute($enst);
        	@rows1 = $sth_db_supp_geneinfo->fetchrow_array();
		$gene_name = $rows1[3]; if (!$gene_name) { $gene_name = ""; }
	        	$gene_desc = $rows1[4]; if (!$gene_desc) { $gene_desc = ""; }
        		$ensfm = $rows1[5]; if (!$ensfm) { $ensfm = ""; }
	        	$fam_desc = $rows1[6]; if (!$fam_desc) { $fam_desc = ""; }
        		$gene_status = $rows1[7]; if (!$gene_status) { $gene_status = ""; }
	        	$fam_size = $rows1[8]; if (!$fam_size) { $fam_size = ""; }
	        	$kaks_mouse = $rows1[9]; if (!$kaks_mouse) { $kaks_mouse = ""; }
		        $kaks_macaque = $rows1[10]; if (!$kaks_macaque) { $kaks_macaque = ""; }
        		$mim_status = $rows1[11]; if (!$mim_status) { $mim_status = ""; }

		unless($original_seq eq ""){
			$comparison_file = "$tmp\/$pid\_$enst\_$hit_id\_comparison.fasta";
		}
		#if indel is in coding region
		if ($region eq "CDS"){
			my @coordInfo = split /,/, $coord;
			my $dis = int($coordInfo[2]) - int($coordInfo[1]);
			my $allelLength = length $coordInfo[4];
			#3n indels, no prediction so far
			#updated on Feb 4, 2013. Now we have predictions for 3n indels
			if ($subst ne "AA_INSERTION" && $subst ne "AA_DELETION" && $subst ne "FRAMESHIFT" && $subst ne "EARLYSTOP" )
			{
				print CLASSIFICATION_FILE "$coord\t$hit_id\t$ensg\t$enst\t$subst\t$coords_change\t$nu_change\t$aa_change\t$region\t$comparison_file\t$indel_span\t$indel_loc_pct\t$visual_tx\t$causes_NMD\t$small_indel_warning\t$repeat_change\t$gene_name\t$gene_desc\t$ensfm\t$fam_desc\t$gene_status\t$fam_size\t$kaks_mouse\t$kaks_macaque\t$mim_status\t$genome1000_annot\t\t\t\n";
			}
			elsif ( $allelLength % 3 ==0 or ( $dis % 3 ==0 and $dis!=0)){
				print CLASSIFICATION_FILE "$coord\t$hit_id\t$ensg\t$enst\t$subst\t$coords_change\t$nu_change\t$aa_change\t$region\t$comparison_file\t$indel_span\t$indel_loc_pct\t$visual_tx\t$causes_NMD\t$small_indel_warning\t$repeat_change\t$gene_name\t$gene_desc\t$ensfm\t$fam_desc\t$gene_status\t$fam_size\t$kaks_mouse\t$kaks_macaque\t$mim_status\t$genome1000_annot\t$pred\t$confidence\t$classificationPath\n";
				}
				#non 3n indels
				else{
				print CLASSIFICATION_FILE "$coord\t$hit_id\t$ensg\t$enst\t$subst\t$coords_change\t$nu_change\t$aa_change\t$region\t$comparison_file\t$indel_span\t$indel_loc_pct\t$visual_tx\t$causes_NMD\t$small_indel_warning\t$repeat_change\t$gene_name\t$gene_desc\t$ensfm\t$fam_desc\t$gene_status\t$fam_size\t$kaks_mouse\t$kaks_macaque\t$mim_status\t$genome1000_annot\t$pred\t$confidence\t$classificationPath\n";
				}
				
			}
		#otherwise no prediction
		else {
			print CLASSIFICATION_FILE "$coord\t$hit_id\t$ensg\t$enst\t$subst\t$coords_change\t$nu_change\t$aa_change\t$region\t$comparison_file\t$indel_span\t$indel_loc_pct\t$visual_tx\t$causes_NMD\t$small_indel_warning\t$repeat_change\t$gene_name\t$gene_desc\t$ensfm\t$fam_desc\t$gene_status\t$fam_size\t$kaks_mouse\t$kaks_macaque\t$mim_status\t$genome1000_annot\t\t\t\n";
			}
		#unlink ("$original_seq");
		#unlink ("$modified_seq");
	}
}
close(CLASSIFICATION_FILE);

#create condensed Classification file (with only one transcript) in case user selects that option.
open( CLASSIFICATION_FILE, "$tmp/$pid.classification_file" )
  || die("Cannot open classification file");
open (CONDENSED_CLASSIFICATION_FILE, ">$tmp/$pid.condensed_classification_file");
my %seen_index = ();
while (<CLASSIFICATION_FILE>){
	chomp;
	my @elts = split /\t/, $_;
	my $coords = $elts[0];
	my $num_elts = num_elts_in_row($_);
	if ($seen_index{$coords} && $seen_index{$coords} ne ""){
		my $prev_row = $seen_index{$coords};
		my $prev_num_elts = num_elts_in_row($prev_row);
		if ($num_elts > $prev_num_elts){
			$seen_index{$coords} = $_;
		}
	}
	else{
		$seen_index{$coords} = $_;
	}
}

foreach my $key (keys %seen_index){
	print CONDENSED_CLASSIFICATION_FILE $seen_index{$key},"\n";
}
close (CONDENSED_CLASSIFICATION_FILE);
close (CLASSIFICATION_FILE);

#read from combined file and print table
my $all_transcripts = defined($opt_m) ? $opt_m:0;

if ($all_transcripts == 1){
	open( CLASSIFICATION_FILE, "$tmp/$pid.condensed_classification_file" )|| die("Cannot open classification file");
}
else{
	open( CLASSIFICATION_FILE, "$tmp/$pid.classification_file" )|| die("Cannot open classification file");
}


$heading_tsv = get_output_heading_tsv($output_options);
open( OUTFILETSV,   ">$tmp/$pid\_predictions.tsv" );
print OUTFILETSV "$heading_tsv\n";

my $warningflag = 0;

while (<CLASSIFICATION_FILE>) {
	chomp;
	@elts = split /\t/, $_;
	$coord =  $elts[0];
	$hit_id = $elts[1];
	$ensg = $elts[2];
	$enst = $elts[3];
	$subst = $elts[4];
	$coords_change = $elts[5];
	$nu_change = $elts[6];
	$aa_change = $elts[7];
	$region = $elts[8];
	$comparison_file = $elts[9];
	$indel_span = $elts[10];
	$indel_loc_pct = $elts[11];
	$visual_tx = $elts[12];
	$causes_NMD = $elts[13];
	$small_indel_warning = $elts[14]; 
	$repeat_change = $elts[15];
	$table_row_html = ""; #	$table_row_html = get_output_row_html($output_options,$_);
	$table_row_tsv = get_output_row_tsv($output_options,$_);
        print OUTFILETSV "$table_row_tsv\n";


=pod
=cut
}



#close(FILE);
close(OUTFILETSV);
close (CLASSIFICATION_FILE);

################################################
# SIMNL<2012-05-29>: directories with underscores causes
# missing fields in results file.
# Temporary fix until we determine why this is so
# Perhaps the new python code does something?
my $final_dir = $tmp;
if ($original_tmp ne $tmp) {
    if (! -e $original_tmp) {
	my $res = system("mkdir 0755 -p $original_tmp");
	print "original_tmp mkdir results: $res -> $original_tmp\n";
    }
    opendir(DIR, $tmp);
	my $f="";
    while($f = readdir(DIR) && $f ne "." && $f ne "..") {
	my $testdir = "$tmp/$f";
	if (-d $testdir) {
	    my $newdir = "$original_tmp/$f";
	    my $r = system("mkdir 0755 -p $newdir");
	    print "$newdir: $r\n";
	    mkdir($newdir, 0755);
	    print "$newdir created\n";
	}
    } #end while
    closedir(DIR);

    my $results = system("cp -R $tmp $original_tmp");
    if ($results == 0) {
	$final_dir = $original_tmp;
	system ("rm -rf $tmp");
	system ("rm -rf $parent");
	system ("rm -rf 0755");
    }


}
################################################

print "predictions in $final_dir/$pid\_predictions.tsv\n";
system ("rm -f $final_dir/$final_dir\_snps_chr*bin*.txt");
system ("rm -f $final_dir/$final_dir\.chrfile");
system ("rm -f $final_dir/$final_dir\.*gff");
system ("rm -f $final_dir/$final_dir\_*.fa");
system ("rm -f $final_dir/$final_dir\_old_new_coords_map.txt");
system ("rm -f $final_dir/$final_dir\_repeat_detect_*file.txt");

=pod 
=cut
my $end = time();
printf("The total running time is %.2f seconds\n", $end - $start);



####gina ADD to check the preprocessor error##
if($prepRet!=0)
{
  print '<BR><FONT COLOR="RED">We found one or more lines in the input data got errors!</FONT>';
  print '<BR>Please see the <A HREF="/www/sift/tmp/'.$pid.'/';
  my $name=`basename $all_chr_file`;
  chomp($name);
  print $name .'.err' ;
  print '">file</A> for detail.';
}
###gina ADD END#####
#Changes of titles made by Jing Hu @ Oct 12, 2011
sub get_output_heading_tsv{
        $oo = $_[0] ;
        @elts = split /,/, $oo;
		#$heading_tsv = "Coordinates\tGene ID\tTranscript ID\tSubstitution Type\tRegion\tAmino acid position change\tIndel Location\tNucleotide change\tAmino acid change\tCauses NMD\tRepeat detected\tTranscript visualization";
        $heading_tsv = "Coordinates\tGene ID\tTranscript ID\tSubstitution Type\tRegion\tAmino acid position change\tEffect\tConfidence Score\tClassification Path\t";
        $heading_tsv .= "Nucleotide change\tAmino acid change\tIndel Location\tCauses NMD\tRepeat detected\tTranscript visualization";
        
        @options = ("Gene Name","Gene Desc","Protein Family ID","Protein Family Desc","Transcript Status","Protein Family Size","Ka/Ks (Mouse)","Ka/Ks (Macaque)","OMIM Disease", "1000 Genomes");


        for ($i = 0 ; $i < scalar @elts; $i++){
                if ($elts[$i] eq "1"){
                        $heading_tsv.= "\t" . $options[$i];
                }
        }
        
        return $heading_tsv;
}


#Changes of html titles made by Jing Hu @ Oct 12, 2011
sub get_output_heading_html{
        $oo = $_[0] ;
        @elts = split /,/, $oo;
        $heading_html = "Coordinates\tGene ID\tTranscript ID\tSubstitution Type\tRegion\t<a href=\"\/www\/indels_help2.html#Amino_acid_position_change\" target=\"_blank\">Amino acid position change</a>";
		$heading_html.= "\tEffect\tConfidence Score\t<a href=\"\/www\/indels_help2.html#classificationPath\" target=\"_blank\">Classification Path</a>\t";
		$heading_html.= "<a href=\"\/www\/indels_help2.html#Nucleotide_change\" target=\"_blank\">Nucleotide change</a>\t<a href=\"\/www\/indels_help2.html#Amino_acid_change\" target=\"_blank\">Amino acid change</a>\t<a href=\"\/www\/indels_help2.html#Protein_sequence_change\" target=\"_blank\">Protein Sequence Change</a>\t<a href=\"\/www\/indels_help2.html#Indel_location\" target=\"_blank\">Indel location</a>\t<a href=\"\/www\/indels_help2.html#NMD\" target=\"_blank\">Causes Nonsense Mediated Decay (NMD)</a>\t<a href=\"\/www\/indels_help2.html#repeat\" target=\"_blank\">Repeat detected</a>";
	
	#$heading_html =
#"Coordinates\tGene ID\tTranscript ID\tSubstitution Type\tRegion\t<a href=\"\/www\/indels_help2.html#Amino_acid_position_change\" target=\"_blank\">Amino acid position change</a>\t<a href=\"\/www\/indels_help2.html#Indel_location\" target=\"_blank\">Indel location</a>\t<a href=\"\/www\/indels_help2.html#Nucleotide_change\" target=\"_blank\">Nucleotide change</a>\t<a href=\"\/www\/indels_help2.html#Amino_acid_change\" target=\"_blank\">Amino acid change</a>\t<a href=\"\/www\/indels_help2.html#Protein_sequence_change\" target=\"_blank\">Protein Sequence Change</a>\t<a href=\"\/www\/indels_help2.html#NMD\" target=\"_blank\">Causes Nonsense Mediated Decay (NMD)</a>\t<a href=\"\/www\/indels_help2.html#repeat\" target=\"_blank\">Repeat detected</a>";
	@options = ("Gene Name","Gene Desc","Protein Family ID","Protein Family Desc","Transcript Status","Protein Family Size","Ka/Ks (Mouse)","Ka/Ks (Macaque)","OMIM Disease", "1000 Genomes");



        for ($i = 0 ; $i < scalar @elts; $i++){
                if ($elts[$i] eq "1"){
			if ($i == 9) {
				$heading_html .= "\t<a href=\"\/www\/indels_help2.html#1000genomes\" target=\"_blank\">$options[$i]</a>";
			} else {
                        	$heading_html.= "\t$options[$i]";
			}
                }
        }
        
        return $heading_html;
}


sub get_output_row_tsv {
        $oo = $_[0];
        $classification_line = $_[1];
	chomp $classification_line;
        @elts = split /\t/, $classification_line;
        $coord =  $elts [0]; if (!$coord) { $coord = ""; } else { chomp($coord); }
        $hit_id = $elts[1];
        $ensg = $elts[2]; if (!$ensg) { $ensg = ""; } else { chomp($ensg); }
        $enst = $elts[3]; if (!$enst) { $enst = ""; } else { chomp($enst); }
        $subst = $elts[4]; if (!$subst) { $subst = ""; } else { chomp($subst); }
        $coords_change = $elts[5]; if (!$coords_change) { $coords_change = ""; } else { chomp($coords_change); }
        $nu_change = $elts[6]; if (!$nu_change) { $nu_change = ""; } else { chomp($nu_change); }
        $aa_change = $elts[7]; if (!$aa_change) { $aa_change = ""; } else { chomp($aa_change); }
        $region = $elts[8]; if (!$region) { $region = ""; } else { chomp($region); }
        $comparison_file = $elts[9];
        $indel_span = $elts[10];
        $indel_loc_pct = $elts[11]; if (!$indel_loc_pct) { $indel_loc_pct = ""; } else { chomp($indel_loc_pct); }
        $visual_tx = $elts[12]; if (!$visual_tx) { $visual_tx = ""; } else { chomp($visual_tx); }
        $causes_NMD = $elts[13]; if (!$causes_NMD) { $causes_NMD = ""; } else { chomp($causes_NMD); }
        $small_indel_warning = $elts[14]; if (!$small_indel_warning) { $small_indel_warning = ""; } else { chomp($small_indel_warning); }


        
	$repeat_change = $elts[15]; if (!$repeat_change) { $repeat_change = ""; } else { chomp($repeat_change); }
		
		#Changes of column orders made by Jing Hu @ Oct 12, 2011
		
		#$table_row_tsv = "$coord $small_indel_warning\t$ensg\t$enst\t$subst\t$region\t$coords_change\t$indel_loc_pct\t$nu_change\t$aa_change\t$causes_NMD\t$repeat_change\t$visual_tx";
		$table_row_tsv = "$coord $small_indel_warning\t$ensg\t$enst\t$subst\t$region\t$coords_change";
	
		#added by Jing Hu, June 13, 2011 to include predictions

    	my $pred = $elts[26];
    	my $confidence = $elts[27];
	my $classificationPath = $elts[28];# if ($subst eq "EARLYSTOP") { print "EARLYSTOP: $pred, $confidence, $classificationPath\n"; }
	if ($pred && $confidence && $classificationPath) { 
		$table_row_tsv.= "\t$pred";
		$table_row_tsv.= "\t$confidence";
		$table_row_tsv.= "\t$classificationPath";
	} else { $table_row_tsv.= "\t\t\t"; }
		$table_row_tsv.="\t$nu_change\t$aa_change\t$indel_loc_pct\t$causes_NMD\t$repeat_change\t$visual_tx";
	
        @elts2 = split /,/, $oo;
        for ($i = 0 ; $i < scalar @elts2; $i++){
	    # Bug fix: 2014.01.13 - chomp returns success (0) or failure, so 
                if ($elts2[$i] eq "1"){
#                if (chomp($elts2[$i]) eq "1"){
			$new_cell = $elts[16+$i];
			if (!$new_cell) { $new_cell = ""; } else { chomp $new_cell; }
                        $table_row_tsv.= "\t$new_cell";
                }
        }
        
        

        return $table_row_tsv;

}

sub get_output_row_html {
        $oo = $_[0];
        $classification_line = $_[1];
        chomp $classification_line;
        @elts = split /\t/, $classification_line;


        $coord =  $elts [0];
        $hit_id = $elts[1];
        $ensg = $elts[2];
        $enst = $elts[3];
        $subst = $elts[4];
        $coords_change = $elts[5];
        $nu_change = $elts[6];
        $aa_change = $elts[7];
        $region = $elts[8];
        $comparison_file = $elts[9];
        $indel_span = $elts[10];
        $indel_loc_pct = $elts[11];
        $visual_tx = $elts[12];
        $causes_NMD = $elts[13];
        $small_indel_warning = $elts[14];
		$repeat_change = $elts[15];
        
        #Changes of column orders made by Jing Hu @ Oct 12, 2011
        
        #$table_row_html ="$coord\t$ensg\t$enst\t$subst\t$region\t$coords_change\t$indel_loc_pct\t$nu_change\t$aa_change\t$comparison_file\t$causes_NMD\t$repeat_change";
        $table_row_html ="$coord\t$ensg\t$enst\t$subst\t$region\t$coords_change";
        #added by Jing Hu, June 13, 2011 to include predictions
        my $pred = $elts[26];
        my $confidence = $elts[27];
		my $classificationPath = $elts[28];
	if ($pred && $confidence && $classificationPath) {
		$table_row_html.= "\t$pred";
		$table_row_html.= "\t$confidence";
		$table_row_html.= "\t$classificationPath"; } else { $table_row_html.= "\t\t\t"; }
        $table_row_html.="\t$nu_change\t$aa_change\t$comparison_file\t$indel_loc_pct\t$causes_NMD\t$repeat_change";
        
	@elts2 = split /,/, $oo;
        for ($i = 0 ; $i < scalar @elts2; $i++){
                if ($elts2[$i] eq "1"){
		    if (not defined $elts[16+$i]) {
			$table_row_html.= "\t&nbsp;";
		    } else {
			my $element = $elts[16+$i];
			$element =~ s/^\s+//;
			$element =~ s/\s+$//;
			if ($element eq "") {
			    $element = "&nbsp;"
			}
			$table_row_html.= "\t$element";
#			$table_row_html.= "\t@elts[16+$i]";
		    }

                }
        }
        
		

        return $table_row_html;

}


sub num_elts_in_row{
	my $row = $_[0];
	my $num_elts;
	chomp $row;
	my @elts = split /\t/, $row;
	foreach my $elt(@elts){
		if ($elt ne ""){
			$num_elts++;
		}
	}	
	return $num_elts;
}


#-------------------------------------------------------------------------
exit(0);

#-------------------------------------------------------------------------
#
# parameter: a string that is the html QUERY_STRING environment
#variable
# returns: an associative array of name/value pairs.  The name is the
#key.
sub parse_query {
	local ($query_string) = @_;
	local ( %ans, @q, $pair );

	#print $query_string;
	# break up into individual name/value lines
	@q = split( /&/, $query_string );

	foreach $pair (@q) {

		# break the name/value pairs up
		# use split rather than regular expressions because the value may
		# have
		#  newlines in it
	    my @split_pair = split( /=/, $pair, 2 );

		# change '+' to ' '
#		$_[1] =~ s/\+/ /g;
	    $split_pair[1] =~ s/\+/ /g;

		# change the escaped characters (has to be after the split on '&'
		# and '=')
#		$_[1] =~ s/%(..)/pack("c",&hextodec("$1"))/eg;
	    $split_pair[1] =~ s/%(..)/pack("c",&hextodec("$1"))/eg;
#		$ans{ $_[0] } = $_[1];
	    $ans{ $split_pair[0] } = $split_pair[1];
	}

	return %ans;
}

#-------------------------------------------------------------------------
# parameter: a hex representation of a number (doesn't need to be a
# string)
# returns: the decimal representation of the number
sub hextodec {
	unpack( "N", pack( "H8", substr( "0" x 8 . shift, -8 ) ) );
}

#-------------------------------------------------------------------------
# $names = &general_parse($ENV{CONTENT_TYPE}, $QUERY_STRING);
# parameters:   CONTENT_TYPE
#               QUERY_STRING
# returns: an associative array of name/value pairs.  The name is the
# key.

# WARNING:  Some of this routine is program-dependent!!!

# CONTENT_TYPE: application/x-www-form-urlencoded
# QUERY_STRING: key1=val1&key2=val2

# CONTENT_TYPE: multipart/form-data; boundary=<boundary>
# QUERY_STRING: <boundary>
#               Content-Disposition: form-data; name="key1"
#               <blank line>
#               val1
#               <boundary>
#               Content-Disposition: form-data; name="key2"
#               <blank line>
#               val2
#               <boundary>

sub general_parse {
	local ( $content_type, $query_string ) = @_;
	local ( %ans, @q, $pair, $loc, $boundary, $temp, $temp1 );

	if ( $content_type eq "application/x-www-form-urlencoded" ) {

		# break up into individual name/value lines
		@q = split( /&/, $query_string );

		foreach $pair (@q) {

			# break the name/value pairs up
			# use split rather than regular expressions because the value
			# may have
			#  newlines in it
		    my @split_pair = split( /=/, $pair, 2 );
		    
			# change '+' to ' '
#			$_[1] =~ s/\+/ /g;
			$split_pair[1] =~ s/\+/ /g;

			# change the escaped characters (must be after the split on '&'
			# and '=')
			$split_pair[1] =~ s/%(..)/pack("c",&hextodec("$1"))/eg;
#			$_[1] =~ s/%(..)/pack("c",&hextodec("$1"))/eg;

#			$ans{ $_[0] } = $_[1];
			$ans{ $split_pair[0] } = $split_pair[1];
		}    #end of foreach $pair

	}    #end of if ($content_type)
	else {
		$loc = index( $content_type, "boundary=" );
		if ( $loc > 0 ) {
			$temp = substr( $content_type, $loc + 9 );

		 #               Why is this necessary? (boundary= doesn't match actual)
			$boundary = "--" . $temp;

			# break up into individual name/value lines
			@q = split( /$boundary/, $query_string );

			foreach $pair (@q) {

				# break the name/value pairs up
				$loc = index( $pair, "name=" );
				$temp = substr( $pair, $loc + 5 );

				#         $loc = index($temp, "\n\n");
				$loc = index( $temp, "\n" );
				$temp1 = substr( $temp, $loc + 2 );

				#   Get rid of stuff after the name; including semicolon if any
				$loc_semi = index( $temp, ";" );
				$loc_eol  = index( $temp, "\n" );
				$loc      = $loc_eol;
				if ( $loc_semi > 0 && $loc_semi < $loc ) {
					$loc = $loc_semi;
				}
				if ( $loc > 0 ) { $temp = substr( $temp, 0, $loc ); }

				#               Get rid of quotes around the name
				$temp =~ s/\"//g;

				#               Still has a trailing whitespace character ...
				$temp =~ s/\s//g;

		  #               Substitute spaces with nothing
		  #               Need to strip leading/ending whitespace off of $temp1,
		  #               but be careful not to strip off internal CRs
		  #               MAC file lines end in just \r, no \n, so makelis won't
		  # find all
		  #               of the sequences; DOS file lines end in \r\n, UNIX in
		  #\n.
		  #               Change \r\n to \n and then \r to \n
#######PROGRAM -SPECIFIC!!!!!!!######################
		 #In my case, I want to keep the newlines in fields which have "file" or
		 # 'seq"
		 # and remove newlines everywhere else.
		 #if ( $temp =~ /file/ || $temp =~ /seq/ || $temp =~ /subst/ ) {
				$temp1 =~ s/\r\n/\n/g;
				$temp1 =~ s/\r/\n/g;

				#}

			 # for other variables that are NOT files or file-like, remove extra
			 #whitespace
			 #else { $temp1 =~ s/\s//g; }
				if ( $temp ne "" ) { $ans{$temp} = $temp1; }
			}    # end of foreach
		}    #end of if loc > 0
		else {
			print "Cannot parse\n";
			print "content_type=$content_type\n";
			print "query_string=$query_string\n";
		}
	}
	return %ans;

	#print "</PRE>";
}    # end of general_parse

# returns hash for a file, 2nd field is the key and the 3rd field
# is the value 4th field, is the delimiter
sub make_hash {
	my ($file) = @_;
	my %hash;
	open( HASH, $file ) || die "can't open $file";
	my $line;
	while ( $line = <HASH> ) {
		chomp($line);
		if ( exists( $hash{$line} ) ) {
			$hash{$line}++;
		}
		else {
			$hash{$line} = 1;
		}
	}
	close(HASH);
	return (%hash);
}
=pod
sub update_IP_logfile {
	my ( $queuefile, $IP_address ) = @_;

	$lockqueuefile = "$queuefile.lock";

	# lockfile will wait until it can lock the file
	`./lockfile $lockqueuefile`;

	# append the address and command to the queue file
	open( FILE, ">>$queuefile" );
	print FILE "$IP_address\n";
	close(FILE);

	chmod( 0664, $queuefile );

	# remove the lock file
	unlink($lockqueuefile);

}
=cut

