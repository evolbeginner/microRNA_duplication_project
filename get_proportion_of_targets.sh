#! /bin/bash

cd $1 > /dev/null

if [ ! -z $2 ]; then
	dup_sin_dir=$2
else
	dup_sin_dir=~/tools/self_bao_cun/search_4_dupli_and_single/TAIR9
fi

#duplicates_list=$dup_sin_dir/duplicates_noMC.list
duplicates_list=$dup_sin_dir/e_values/duplicates.e-30_noMC.list
singletons_list=$dup_sin_dir/singletons_noMC.list

##############################################################
num_of_dup=`wc -l $duplicates_list | grep -o '^[0-9]\+'`
num_of_sin=`wc -l $singletons_list | grep -o '^[0-9]\+'`

for i in target.{psRNAtarget,tapir,UEA,2_outof_3}.list; do
	unset num_of_targets
	declare -a num_of_targets
	echo $i
	num_of_targets_in_dup=`awk '{print $1}' $i  | xargs -i grep {} $duplicates_list | wc -l`;
	num_of_targets_in_sin=`awk '{print $1}' $i  | xargs -i grep {} $singletons_list | wc -l`;
	num_of_total_targets=`wc -l $i | grep -o '^[0-9]\+'`
	for j in dup sin; do
		eval num1=\${num_of_targets_in_$j}
		eval num_of_type=\${num_of_$j}
		num_of_targets=(${num_of_targets[@]} $num1)
		echo $num1 $num_of_type $num_of_total_targets
		echo "scale=3; $num1/$num_of_total_targets" | bc
	done
	chitest.pl ${num_of_targets[0]} $num_of_dup ${num_of_targets[1]} $num_of_sin
	echo
done

cd - >/dev/null

