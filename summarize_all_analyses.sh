#! /bin/bash

dir=$1
list_arr=($dir/either.list $dir/both.list)
expression_result_dir="expression/expr_result"

################################################################################
for list in ${list_arr[@]}; do
	name1=`basename \`dirname $list\``
	name2=`basename $list`
	name2=${name2%.*list}
	echo $name1, $name2, $list
	#array=(¡°${array[@]}¡± $new_element)
	output=("${output[@]}" "expression/expr_result/$name1.$name2.expr")
	perl expression.pl --pair_list_file $list --pair_seperator='-' | awk '{print $3}' > $expression_result_dir/$name1.$name2.expr
done


ttest.pl ${output[0]} ${output[1]}
echo
wilcox_test.pl ${output[0]} ${output[1]}
echo

perl ka_PD_VS_nonPD.pl -in1 ${list_arr[0]} -in2 ${list_arr[1]} -ka_file ../official/kaks/ATH_WGD.kaks -ka_file ../official/kaks/ATH_TD_mayer.kaks -ka_file ../official/kaks/ATH_SSD_my.kaks -kaks_swi -seq_list_file_seperator '-' --is_split_tab_first
echo


