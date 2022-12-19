#!/bin/bash
####################################################################
# Newborn Screening (NBS) Pipeline
#
# Texas Department of State Health Services
#
# Samantha Marcellus, Kyle Long, Rachel Lee
#
# For questions, please contact Kyle Long (kyle.long@dshs.texas.gov
# or Samantha Marcellus (samantha.marcellus@dshs.texas.gov)
###################################################################

#gather user input here
locus="$1"
gene="$2"
var="$3"
out_name="$4"
transcript="$5"
np_id="$6"
date_dir="$7"
run_dir="$8"
user="$9"
date="${10}"

#Set base dir
base_dir=/work/NBS/

#start pipeline
printf "\nGene: $gene"
printf "\nLocus: $locus"
printf "\nVariant: $var"

#python3.7 /home/klong/bb_django/bluebonnet/start_query_tmp.py $gene $var

#Specify analysis folder
prefix=$run_dir/$out_name
#latest_dir=$base_dir/outputs/latest/
sift_name=$(echo $out_name | sed "s/_/-/g")

#ensembl
printf "\nGathering information from ensembl...\n"
ensembl_url="https://rest.ensembl.org/variant_recoder/human/"$locus":"$var"?"
if wget -q -O /dev/null --header='Content-type:application/json' "$ensembl_url"; then
	wget -q --header='Content-type:application/json' "$ensembl_url"  -O $prefix\_ensembl.out
	rsid_backup=$(cat $prefix\_ensembl.out | jq -r '.[0] | to_entries[] | .value["id"][0]')
	norm_var=$(cat $prefix\_ensembl.out | jq -r '.[0] | to_entries[] | .value["hgvsc"]'| grep "$locus" | cut -d ":" -f2 | cut -d '"' -f1 )
else
	rsid_backup="null"
	norm_var="null"
fi

#check for difference in normalized variant and input variant
printf "\nChecking for normalized variant format..."
if [[ "$norm_var" = "$var" || "$norm_var" = "null" || "$norm_var" = "" ]] ; then
        printf "\nVariant format is unchanged\n"
	orig_var=$var
else
	printf "\nVariant format has been normalized\n"
	orig_var=$var
	var=$norm_var
	printf "Input: $orig_var\nNormalized: $var \n"
fi

#dbSNP
echo "$locus":"$var" > $run_dir/var_serv.tmp
test_rsid=$(python3.6 /work/NBS/dbSNP/spdi_batch.py -i $run_dir/var_serv.tmp -t HGVS_RS | tail -1 | cut -f3)
if [[ "$test_rsid" = "no rs found" || "$test_rsid" == *"ERROR"* ]] ; then
	rsid=$rsid_backup
else
	rsid="rs"$test_rsid
fi

if [ "$rsid" = "null" ] ; then
        printf "\nThis mutation does not have an entry in dbSNP\n"
        touch $prefix\_dbSNP.txt
        rsid="NA"
        dbSNP_status=0
else
        printf "\nGathering information from dbSNP...\n"
        dbsnp_url="https://rest.ensembl.org/variation/human/"$rsid"?"
        wget -q --header='Content-type:application/json' "$dbsnp_url"  -O - | jq . > $prefix\_dbSNP.txt
        Rscript $base_dir/dbSNP/dbSNP.R $rsid >> $prefix\_dbSNP.txt
	esearch -db snp -query "$rsid" | efetch -format json | jq -r '.primary_snapshot_data.allele_annotations | .. | .frequency? | select( . != null ) | .[] | select(.observation | .inserted_sequence!=.deleted_sequence) | .study_name, .study_version, .observation.deleted_sequence, .observation.inserted_sequence, .allele_count, .total_count' | paste - - - - - - > $prefix\_dbSNP_freqs.txt
        dbSNP_status=1
fi

#GnomAD (and clinvar)
printf "\nGathering information from gnomad...\n" #This step gathers some information from ClinVar as well
python3.6 $base_dir/gnomad_python_api/gnomad_api_cli.py -filter_by=gene_symbol -search_term=$gene -dataset="gnomad_r3" -reference_genome='GRCh38' > $run_dir/gnomad_r3_raw.tmp
gnomad_check=$(wc -l $run_dir/gnomad_r3_raw.tmp | cut -d ' ' -f1)
if [ "$gnomad_check" = 0 ]; then
	printf "\nWARNING: gnomad did not return any information on this gene. An error likely occurred when communicating with the database. This query should be run again."
	echo "Mutation not found" > $prefix\_gnomad.txt
	echo "Mutation not found" > $run_dir/clinvar.tmp
else
	python3.6 $base_dir/gnomad_python_api/gnomad_parser.py $run_dir/gnomad_r3_raw.tmp $var > $prefix\_gnomad.txt
	python3.6 $base_dir/gnomad_python_api/clinvar_parser.py $run_dir/gnomad_r3_raw.tmp $var > $run_dir/clinvar.tmp
fi

#Exac
printf "\nGathering information from ExAC..."
printf "\nNOTE: ExAC is only availble for GRCh37, so results will be for that reference only\n"
python3.6 $base_dir/gnomad_python_api/exac_api_cli.py -filter_by=gene_symbol -search_term=$gene -dataset="exac" -reference_genome='GRCh37' > $run_dir/exac_raw.tmp
exac_check=$(wc -l $run_dir/exac_raw.tmp | cut -d ' ' -f1)
if [ "$exac_check" = 0 ]; then
        printf "\nWARNING: exac did not return any information on this gene. An error likely occurred when communicating with the database. This query should be run again."
        echo "Mutation not found" > $prefix\_exac.txt
else
	python3.6 $base_dir/gnomad_python_api/gnomad_parser.py $run_dir/exac_raw.tmp $orig_var > $prefix\_exac.txt 
fi

#Collect status of outputs
if grep -q "Mutation not found" $prefix\_gnomad.txt; then gnomad_status=0; else gnomad_status=1; fi
if grep -q "Mutation not found" $run_dir/clinvar.tmp; then clinvar_status=0; else clinvar_status=1; fi

#Check if the mutation was found in any of the first 3 DBs
if [ "$gnomad_status" = 1 ] ; then
	info_file=$prefix\_gnomad.txt
elif [ "$clinvar_status" = 1 ] ; then
	info_file=$run_dir/clinvar.tmp
else
	info_file=0
fi

#fail and exit
if [[ "$info_file" = 0 && "$dbSNP_status" = 0 ]]; then
	printf "\nThe lack of information in these databases hinders downstream analysis. There will be no results for this mutation search\n"
	cp $prefix.log $date_dir 
	#cp $prefix.log $latest_dir
	echo "First fail"
	python3.7 $base_dir/excel_parse/excel_fill_fail.py $prefix $gene $var $locus $user $date
	cp $prefix\_variant_worksheet.xlsx $date_dir
	#cp $prefix\_variant_worksheet.xlsx $latest_dir
	exit
fi

#Gather nucleotide mutation from clinvar/gnomad output
if [ "$dbSNP_status" = 1 ] ; then
	chr=$(grep "seq_region_name" $prefix\_dbSNP.txt | sed 's/ //g' | cut -d ":" -f2 | sed 's/.$//' | sed 's/"//g')
        pos=$(grep '"start":' $prefix\_dbSNP.txt | sed 's/ //g' | cut -d ":" -f2 | sed 's/.$//')
        ref_al=$(grep "allele_string" $prefix\_dbSNP.txt | sed 's/ //g' | cut -d ":" -f2 | sed 's/.$//' | sed 's/"//g' | cut -d "/" -f1)
        alt_al=$(grep "allele_string" $prefix\_dbSNP.txt | sed 's/ //g' | cut -d ":" -f2 | sed 's/.$//' | sed 's/"//g' | cut -d "/" -f2)
	#hgvsp=$(grep "$np_id" $prefix\_dbSNP.txt |  cut -d ':' -f 2 |  sed 's/,//' | sed 's/"//g')
else
	chr=$(grep "variant_id" $info_file | cut -d '-' -f 1 | cut -d ' ' -f 2 )
	pos=$(grep "variant_id" $info_file | cut -d '-' -f 2 )
	ref_al=$(grep "variant_id" $info_file | cut -d '-' -f 3 )
	alt_al=$(grep "variant_id" $info_file | cut -d '-' -f 4 | cut -d ',' -f1 )
	#hgvsp=$(grep "hgvsp:" $info_file | cut -d ' ' -f 2)
fi

#hgvsp
hgvsp=$(mutalyzer_name_checker "$locus:$var" | jq -r '.protein | .description' | cut -d ":" -f 2 | tail -1 | sed 's/(//g' | sed 's/)//g')


#get lengths
len_ref=${#ref_al}
len_alt=${#alt_al}
if [ $len_ref -gt $len_alt ] ; then #deletion
	end=$((pos+len_ref-1))
	ucsc_end=$((pos+len_ref))
else #insertion or SNP
	end=$pos
	ucsc_end=$((pos+1))
fi

printf "\nInformation for downstream tools..."
printf "\nChromosome: $chr"
printf "\nPosition: $pos"
printf "\nReference Allele: $ref_al"
printf "\nAlternate Allele: $alt_al"
#printf "\nrsID: $rsid\n"

#fail and exit
if [ "$chr" = '' ] || [ "$ref_al" = '' ] || [ "$alt_al" = '' ] || [ "$pos" = '' ]; then
#if [[ "$info_file" = 0 && "$dbSNP_status" = 0 ]]; then
	echo "Second fail"
        printf "\nThe lack of information in these databases hinders downstream analysis. There will be no results for this mutation search\n"
	cp $prefix.log $date_dir
        #cp $prefix.log $latest_dir
        python3.7 $base_dir/excel_parse/excel_fill_fail.py $prefix $gene $var $locus $user $date
        cp $prefix\_variant_worksheet.xlsx $date_dir
        #cp $prefix\_variant_worksheet.xlsx $latest_dir
	exit
fi


#PolyPhen2
printf "\nGetting PolyPhen2 scores...\n"
cd $base_dir/ensembl-vep
./vep --id "$locus:$var" -o $prefix\_polyphen2.txt --database --polyphen b --no_stats
cd $base_dir

#ClinVar
#This second clinvar search is needed, because the clinvar search through gnomad is limited
if [ "$clinvar_status" = 0 ] ; then
	printf "\nThis mutation does not have an entry in ClinVar\n"
	echo "N/A" > $prefix\_clinvar.txt 
else
	printf "\nGathering information from ClinVar...\n"
	esearch -db clinvar -query "$locus:$var" | efetch -format variationid |  xtract -pattern ClinicalAssertion -tab " | " -def "N/A" -element "ClinicalAssertion@ID", Description, ReviewStatus, "ClinicalAssertion@DateCreated", "ClinicalAssertion@DateLastUpdated", "ClinicalAssertion@DateSubmitted", "Interpretation@DateLastEvaluated", Comment > $prefix\_clinvar.txt

fi

#Splice Sites
printf "\nCopying saved splice site information...\n"
cp $base_dir/saved_data/$gene\_fsplice.txt $prefix\_fsplice.txt

#UCSC
printf "\nGathering information from UCSC (phyloP)...\n"
UCSC_string='https://api.genome.ucsc.edu/getData/track?genome=hg38;chrom=chr'$chr';start='$pos';end='$ucsc_end';track=phyloP100way'
curl -s -L $UCSC_string | python -c "exec(\"import sys, json\nout=list(json.load(sys.stdin)['chr$chr'])\nfor i in out:\n\tprint(i)\")" > $prefix\_UCSC.txt

#CADD
printf "\nRunning CADD...\n"
printf "$chr\t$pos\t.\t$ref_al\t$alt_al" > $run_dir/CADD_vcf_tmp.vcf
source /home/klong/miniconda3/etc/profile.d/conda.sh
conda activate CADD
$base_dir/CADD-scripts/CADD.sh -g GRCh38 -o $prefix\_CADD.txt.gz $run_dir/CADD_vcf_tmp.vcf
conda deactivate
gunzip $prefix\_CADD.txt.gz

#SIFT 
printf "\nRunning SIFT...\n"
if [ $len_ref = $len_alt ] ;
then
	cp $base_dir/SIFT/vcf_template.vcf $prefix\_SIFT.vcf
	printf "$chr\t$pos\t$rsid\t$ref_al\t$alt_al\t30\tPASS\tNS=1" >> $prefix\_SIFT.vcf
	java -jar $base_dir/SIFT/SIFT4G_Annotator.jar -c -i $prefix\_SIFT.vcf -d $base_dir/SIFT/GRCh38.83.chr -r $run_dir -t
else
	if [ $len_ref -gt $len_alt ] ; #deletion
	then
		printf "$chr,$pos,$end,1,$ref_al/$alt_al" > $run_dir/SIFTINDEL_input_tmp.vcf
	else #insertion
		printf "$chr,$pos,$pos,1,$ref_al/$alt_al" > $run_dir/SIFTINDEL_input_tmp.vcf
	fi
	$base_dir/SIFT/SIFTINDEL/bin/SIFT_exome_indels.pl -i $run_dir/SIFTINDEL_input_tmp.vcf -c hs38 -d $base_dir/SIFT/SIFT_INDEL_HG38 -o $run_dir -p $sift_name\-SIFT
fi

#InterVar
#printf "\nRunning InterVar...\n"
#printf "$chr\t$pos\t$end\t$ref_al\t$alt_al" > $output_dir/intervar_input.tmp
#cd $base_dir/InterVar
#python3.6 Intervar.py  -b hg38 -i $output_dir/intervar_input.tmp  --input_type=AVinput  -o $prefix\_intervar
#cd $base_dir

#PROVEAN
printf "\nRunning Provean...\n"
#printf "$hgvsp"
if [ "$hgvsp" = 'null' ] || [ "$hgvsp" = "" ] || [ "$hgvsp" = 'p.=' ]; then
	printf "Provean evaluates protein mutations, and this mutation does not affect the CDS."
	echo "N/A: Mutation does not affect the protein." > $prefix\_provean.txt
	hgvsp="N/A"
elif [[ "$hgvsp" == *"fs"* ]]; then
	printf "PROVEAN does not evaluate frameshift mutations."
	echo "N/A: PROVEAN does not evaluate frameshift mutations." > $prefix\_provean.txt
else
	python3.6 $base_dir/hgvs/translate.py "$locus:$hgvsp" | cut -d ":" -f 2 | cut -d '.' -f 2 > $run_dir/provean_input.tmp
	provean.sh -q $base_dir/saved_data/$gene.uniprot.fasta --supporting_set $base_dir/saved_data/$gene.sss -v $run_dir/provean_input.tmp --num_threads 16 -V > $prefix\_provean.txt
fi

#clean up
printf "\nCleaning output directory...\n"
mkdir $run_dir/temp_dir
mv $run_dir/*tmp* $run_dir/temp_dir/
rm -f $prefix\_SIFT.vcf
#[[ -f $prefix\_intervar.hg38_multianno.txt ]] && mv $prefix\_intervar.hg38_multianno.txt $prefix\_intervar.hg38_multianno.txt.full
#[[ -f $prefix\_intervar.hg38_multianno.txt.intervar ]] && mv $prefix\_intervar.hg38_multianno.txt.intervar $prefix\_intervar.txt
echo $sift_name
[[ -f $prefix\_SIFT_SIFTpredictions.vcf ]] && mv $prefix\_SIFT_SIFTpredictions.vcf $prefix\_SIFT.txt
[[ -f $run_dir/$sift_name\-SIFT/$sift_name\-SIFT_predictions.tsv ]] && cp $run_dir/$sift_name\-SIFT/$sift_name\-SIFT_predictions.tsv $prefix\_SIFTindel.txt
mkdir $run_dir/other_outputs
#mv $run_dir/*hg38_multianno* $output_dir/other_outputs
[[ -d $run_dir/$sift_name\-SIFT ]] && mv $run_dir/$sift_name\-SIFT $run_dir/other_outputs
[[ -f $prefix\_SIFT_SIFTannotations.xls ]] && mv $prefix\_SIFT_SIFTannotations.xls $run_dir/other_outputs

#Parse raw outputs to excel file
printf "\nParsing information and creating Excel worksheet...\n"
for i in $run_dir/*.txt; do
	echo $i >> $prefix.all;
	cat $i >> $prefix.all;
	printf "\n************************************************************************************************************************\n\n\n" >> $prefix.all;
done;

python3.6 $base_dir/excel_parse/excel_fill.py $prefix $rsid $gene $var $locus $pos $end $chr $ref_al $alt_al $transcript $hgvsp $user $date
if [ -f $prefix\_variant_worksheet.xlsx ]; then
	#cp $prefix\_variant_worksheet.xlsx $date_dir
	#cp $prefix\_variant_worksheet.xlsx $latest_dir
	cp $prefix.log $date_dir
        #cp $prefix.log $latest_dir
else
	cp $prefix.log $date_dir
	#cp $prefix.log $latest_dir
fi

