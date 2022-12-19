#!/usr/local/bin/perl -w

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.
use strict;
my $chrfile = $ARGV[0];

my $tempfile = $chrfile . ".temp";
my $orgfile = $chrfile . ".temp2";
# Pauline March 28,2011
# change pid to just chrfile_name

open (TEMPFILE, ">$tempfile") || die ("cannot open tmp file $tempfile");
open (ORGFILE, ">$orgfile") || die ("cannot open tmp file $orgfile");
while (<>){
	chomp;
	unless (/APPLICATION/) {
	$_ =~ s/^\s+|\s+$//g;
	my $chr = "";
	my $start = "";
	my $stop = "";
	my $orn = "";
	my $allele = "";
	my $comment = "";
	my @elts = split /,/, $_;
	$chr = $elts[0];
	$start = $elts[1];
	$stop = $elts[2];
	$orn = $elts[3];
	$allele = $elts[4];
	$comment = $elts[5];
	next if ($chr eq "");
	#rare bad case if stop < start
	if ($start > $stop){
		next;
	}	

	#if insertion
	if ($start == $stop){
		if ($allele =~ /(\w+)/){	#if already in snp_classifier format (no slash): ATGGC
			$allele = uc ($1);
		}
		elsif($allele =~ /^\-*\/(\w+)/){  #if of the form -/ATTGCA -> ATTGCA
			$allele = uc ($1)
		}
		elsif ($allele =~ /^(\w+)\/\-*/){  #if of the form ATTG/- -> ATTG
			$allele = uc ($1);
		}
	}

	#if deletion
	if ($stop - $start >= 1){
		$allele = "-/";
	}
	if ($comment && $comment !~ /\w|\d/){
		print TEMPFILE "$chr,$start,$stop,$orn,$allele\n";
		if ($start == $stop) {
			print ORGFILE "$chr,$start,$stop,$orn,$allele\n";
			}
		else{
			print ORGFILE "$chr,$start,$stop,$orn,/\n";
			}
	} elsif (!$comment) {
		print TEMPFILE "$chr,$start,$stop,$orn,$allele\n";
		if ($start == $stop) {
			print ORGFILE "$chr,$start,$stop,$orn,$allele\n";
			}
		else{
			print ORGFILE "$chr,$start,$stop,$orn,/\n";
			}		
	}
	else{	    
		print TEMPFILE "$chr,$start,$stop,$orn,$allele,$comment\n";
		if ($start == $stop) {
			print ORGFILE "$chr,$start,$stop,$orn,$allele,$comment\n";
			}
		else{
			print ORGFILE "$chr,$start,$stop,$orn,/,$comment\n";
			}
		
	}
	} # end match it is an CONTENT-TYPE APPLICATION line
}
system ("cp $orgfile $chrfile.orginal");
system ("mv $tempfile $chrfile");
