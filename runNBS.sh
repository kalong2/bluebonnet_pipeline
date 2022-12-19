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

#Set base dir
base_dir=/work/NBS/
#user=$(id -u -n)
#user=$(logname)
out_name=$(basename $1 | cut -d '.' -f 1)
#excel_dir=$base_dir/outputs/Variant_Worksheets/$out_name
date_str=$2
#latest_dir=$base_dir/outputs/latest/
date_dir=$base_dir/outputs/$date_str
run_dir=$base_dir/outputs/$date_str/$out_name
echo $run_dir
rm -rf $run_dir
mkdir -p $run_dir
user=$3

declare -A locus_dict=( ["ABCD1"]="NM_000033.4" ["HBB"]="NM_000518.5" ["ACADVL"]="NM_000018.4")
declare -A transcript_dict=( ["ABCD1"]="ENST00000218104" ["HBB"]="ENST00000335295" ["ACADVL"]="ENST00000356839")
declare -A np_dict=( ["ABCD1"]="NP_000024.2" ["HBB"]="NP_000509.1" ["ACADVL"]="NP_000009.1")

while read line; do
	line=$(echo $line | xargs)
	if [ -z "$line" ] ; then
		continue
	fi
	gene=$(echo $line | cut -d ' ' -f1)
	var=$(echo $line | cut -d ' ' -f2)
	locus=${locus_dict[$gene]}
	transcript=${transcript_dict[$gene]}
	np_id=${np_dict[$gene]}

	#Check Variant List
	#check="$gene $var"
	#if grep -q "$check" $base_dir/variant_list.txt ;
	#then
        #	id=$(grep "$check" $base_dir/variant_list.txt | tr -s ' '| cut -d ' ' -f1)
	#
	#else
        #	id=$(wc -l < $base_dir/variant_list.txt)
        #	printf "$id $gene $var\n" >> $base_dir/variant_list.txt
	#fi



	#nohup $base_dir/NBSpipeline.sh $locus $gene $var $id $user > $output_dir/$gene\_$id.log &
	echo "nohup $base_dir/NBSpipeline.sh '$locus' '$gene' '$var' '$out_name' '$transcript' '$np_id' '$date_dir' '$run_dir' '$user' '$date_str'> $run_dir/$out_name.log 2>&1 &" | batch
	python3.7 /home/klong/bb_django/bluebonnet/start_query_tmp.py $gene $var
done < $1

