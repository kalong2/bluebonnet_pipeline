#!/usr/bin/perl -w 

use strict;
#use File::Copy;
use File::Path; # for rmtree

my $int_bin = "./int_bin";
my $fin_bin = "./bin";

my $user_config = "./config_env.txt";
my %user_config_hash = ();

my $regexpfile = "./regexpfile.txt";
my @file_regexp_replace_key = ();
my $placeholder = "PLACEHOLDER";

###### MAIN #######

&make_a_directory($fin_bin, $int_bin); # create the intermediate directory

&read_user_config(); # read in user's configuration

# Create the tmp directory if it does not exist
my $tmp_directory = $user_config_hash{TMP_DIR};
if (defined $tmp_directory && !(-e $tmp_directory)) {
    print "Specified tmp directory does not exist, creating directory\n";
    my $result = system("mkdir $tmp_directory");
    if ($result != 0) {
	die "Creating tmp directory: Unable to create $tmp_directory\n";
    } else {
	$result = system("chmod 775 $tmp_directory");
    }
}



&read_regexp_file($regexpfile); # read file:regexp:replacement:key 
&perform_replacement();

&changePermissions();

my $result = system("rm -rf $int_bin");

print "Environment setup completed\n";

###### SUB-ROUTINES ######
sub changePermissions() {

    my $results = system("chmod 755 ./bin/Classify_SNPs_Indels_jhu_March2013.pl");
    $results = system("chmod 755 ./bin/model_transcript.pl");
    $results = system("chmod 755 ./bin/SIFT_exome_indels.pl");
	$results = system("chmod 755 ./bin/map_coords_to_bin_indels.pl");
	$results = system("chmod 755 ./bin/indel-bin/*.py");
}


sub perform_replacement() {
    for my $a (@file_regexp_replace_key) {
	my ($filename, $regexp, $replace, $key) = split(/:/,$a);

	my $contents = "";

	open (INFILE, "<$int_bin/$filename") or die "Unable to open $int_bin/$filename\n";
	while(my $line = <INFILE>) {
	    chomp $line;
	    if ($line =~ /$regexp/) {
		my $value = $user_config_hash{$key};
		$replace =~ s/$placeholder/$value/;
		$line = $replace;
	    } 
	    $contents = "$contents$line\n";
	}
	close(INFILE);

	open(OUTFILE, ">$fin_bin/$filename") or die "Perform_replacement: Unable to open $fin_bin/$filename\n";
	print OUTFILE $contents;
	close(OUTFILE);

	# Write back to intermediate file because the file may have multiple replacements.
	open(OUTFILE, ">$int_bin/$filename") or die "Perform_replacement: Unable to open $int_bin/$filename\n";
	print OUTFILE $contents;
	close(OUTFILE);
    }
}#end perform_replacement

sub read_regexp_file($) {
    my ($file) = @_; 
    open(FILE, "<$file") or die "Unable to open file $file\n";
    while (my $line = <FILE>) {
	chomp $line;
	if (not ($line =~ /^\s*$/) && not ($line =~ /^\#/)) {
	    $line = trimstr($line);
	    push @file_regexp_replace_key, $line;
	} 
    }
    close(FILE);
}#end read_regexp_file

sub read_user_config() {
    open(USERCONFIG, "<$user_config") or die "readUserConfig: Unable to open $user_config\n";
    while (my $line = <USERCONFIG>) {
	chomp $line;
	my ($key, $value) = split(/=/, $line);
	$key = trimstr($key);
	$value = trimstr($value);
	$user_config_hash{$key} = $value;
    }
    close(USERCONFIG);
}#end read_user_config

sub make_a_directory($$) {
    my ($origdir, $newdir) = @_;
    if (-e $newdir) {
	rmtree($newdir);
    }
    my $result = system("cp -R $origdir $newdir");
    if ($result != 0) {
	die "make_a_directory: Unable to create $newdir\n";
    }
}#end make_a_directory

sub trimstr($) {
    my ($str) = @_;
    $str =~ s/^\s+//;
    $str =~ s/\s+$//;
    return $str;
}#end trimstr

__END__
