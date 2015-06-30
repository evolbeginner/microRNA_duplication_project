#! /bin/bash

declare -A num_of_WGD num_of_SSD num_of_TD
diver_types=(both either neither)
duplicate_types=(WGD SSD TD)
target_dir=$1

##############################################################################
cd $target_dir

for duplicate_type in ${duplicate_types[@]}; do	
	echo $duplicate_type
	for diver_type in ${diver_types[@]}; do
		#eval num_of_$duplicate_type[$diver_type]=`wc -l ${diver_type}_${duplicate_type}.list | grep -o '^[0-9]\+'`
		eval num_of_$diver_type=`wc -l ${diver_type}_${duplicate_type}.list | grep -o '^[0-9]\+'`
	done
	#eval tmp_array=\${num_of_$duplicate_type[@]}
	#for i in ${tmp_array}; do
	#	tmp_array2=(${tmp_array2[@]} $i)
	#done
	num_of_with_microRNA_site=`expr $num_of_either + $num_of_both`
	num_of_all=`expr $num_of_either + $num_of_both + $num_of_neither`
	echo -e "either: $num_of_either\t$num_of_with_microRNA_site"
	echo "scale=3; $num_of_either / $num_of_with_microRNA_site" | bc
	echo -e "with_microRNA_site: $num_of_microRNA_site\t$num_of_all"
	echo "scale=3; $num_of_with_microRNA_site / $num_of_all" | bc
done


