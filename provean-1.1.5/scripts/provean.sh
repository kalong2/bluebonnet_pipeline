#!/bin/bash
#
# Copyright 2012-2014 J. Craig Venter Institute
# This file is part of PROVEAN.  PROVEAN is free software: you may
# redistribute it and/or modify it under the terms of the GNU General Public
# License version 3.  PROVEAN is distributed in the hope that it will be
# useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.  You should have received a copy of the GNU
# General Public License along with this program.  If not,
# see http://www.gnu.org/licenses/gpl.txt.
#
# Name		: provean.sh
# Author	: Yongwook Choi 
# 
##############################################

####################
# CONFIGURATION
####################
# Specify the path to database and program
#
BLAST_DB="/work/NBS/databases/nr_v4/nr"
PSIBLAST="/work/NBS/tools/ncbi-blast-2.4.0+/bin/psiblast"
CD_HIT="/work/NBS/tools/cdhit/cd-hit"
BLASTDBCMD="/work/NBS/tools/ncbi-blast-2.4.0+/bin/blastdbcmd"
# END CONFIGURATION
####################



shopt -s -o nounset
SCRIPT="provean.sh"
SCRIPT_DIR=$(readlink -f $0)
SCRIPT_DIR=${SCRIPT_DIR%/*}

if [ -z "$BLAST_DB" ] ; then
	echo "error: BLAST database name is missing. Please edit provean.sh file to add the name."
	exit 1;
fi

if [ -z "$PSIBLAST" ] ; then
	echo "error: psiblast path is missing. Please edit provean.sh file to add the path."
	exit 1
fi

if [ -z "$CD_HIT" ] ; then
	echo "error: cd-hit path is missing. Please edit provean.sh file to add the path."
	exit 1
fi

if [ -z "$BLASTDBCMD" ] ; then
	echo "error: blastdbcmd path is missing. Please edit provean.sh file to add the path."
	exit 1
fi

QUERY=
VARIATION=
QUIET="--quiet"
SSS=
SAVE_SSS=
VERBOSE=
NUM_THREADS=
TMP_DIR=

# check getopt mode
getopt -T
if [ $? -ne 4 ] ; then 
	echo "error: Requires enhanced getopt, obtain new version."
	exit 1;
fi

OPTSTRING="q:v:Vh"
LOPTSTRING="query:,variation:,save_supporting_set:,supporting_set:,num_threads:,tmp_dir:,verbose,help"
USAGE="PROVEAN v1.1.5

USAGE:
  provean.sh [Options]

Example:
 # Given a query sequence in aaa.fasta file, 
 # compute scores for variations in bbb.var file 
 provean.sh -q aaa.fasta -v bbb.var

Required arguments:
 -q <string>, --query <string>
   Query protein sequence filename in fasta format
 -v <string>, --variation <string>
   Variation filename containing a list of variations:
     one entry per line in HGVS notation,
     e.g.: G105C, F508del, Q49dup, Q49_P50insC, Q49_R52delinsLI

Optional arguments:
 --save_supporting_set <string>
   Saves supporting sequence set infomation into a given filename
 --supporting_set <string>
   Supporting sequence set filename saved with '--save_supporting_set' option above
   (This will save time for BLAST search and clustering.)
 --tmp_dir <string>
   Temporary directory used to store temporary files
 --num_threads <integer>
   Number of threads (CPUs) to use in BLAST search
 -V, --verbose
   Verbosely shows the information about procedure
 -h, --help
   Gives this help message
"

RESULT=$(getopt -n "$SCRIPT" -o "$OPTSTRING" -l "$LOPTSTRING" -- "$@")
if [ $? -ne 0 ] ; then
	# parsing error, show usage
	echo "$USAGE" 
	exit 1
fi

eval set -- "$RESULT"
while [ true ] ; do
	case "$1" in
		-q|--query) 
			shift 
			QUERY="$1"
		;;
		-v|--variation)
			shift
			VARIATION="$1"
		;;
		-V|--verbose)
			QUIET=""
		;;
		--supporting_set)
			shift
			SSS="$1"
		;;
		--save_supporting_set)
			shift
			SAVE_SSS="$1"
		;;
		--tmp_dir)
			shift
			TMP_DIR="$1"
		;;
		--num_threads)
			shift
			NUM_THREADS="$1"
		;;
		-h|--help)
			echo "$USAGE"
			exit 0
		;;
		--)
			shift
			break
		;;
	esac
	shift
done

if [ -z "$QUERY" ] ; then
	echo "error: need query sequence filename" 
	exit 1
fi

if [ -z "$VARIATION" ] ; then
	echo "error: need variation filename"
	exit 1
fi

COMMAND="$SCRIPT_DIR/provean -q $QUERY -v $VARIATION -d $BLAST_DB --psiblast $PSIBLAST --cdhit $CD_HIT --blastdbcmd $BLASTDBCMD $QUIET"

if [ -n "$SAVE_SSS" ] ; then
	COMMAND="$COMMAND --save_supporting_set $SAVE_SSS"
fi

if [ -n "$SSS" ] ; then
	COMMAND="$COMMAND --supporting_set $SSS"
fi

if [ -n "$TMP_DIR" ]; then
	COMMAND="$COMMAND --tmp_dir $TMP_DIR"
fi

if [ -n "$NUM_THREADS" ]; then
	COMMAND="$COMMAND --num_threads $NUM_THREADS"
fi

# run command
$COMMAND

STATUS=$?

exit $STATUS

