#! /bin/bash

e_value=1e-10

############################################################################################################
getopt_cmd="getopt -o 'h' --long help,blast_result_dir:,seq:,seq_file:,search_4_dupli_and_single_prog:,outdir:,e_value: -n 'bmtagger' -- "$@""
PARAM=$($getopt_cmd)
eval set -- "$PARAM"
while true ; do
	case "$1" in
		-h|--help)		show_help;		shift ;;
		--blast_result_dir)	blast_result_dir=$2;	shift 2 ;;
		--search_4_dupli_and_single_prog)
					search_4_dupli_and_single_prog=$2;	shift 2 ;;
		--seq_file|--seq)	seq=$2;			shift 2 ;;
		--e_value)		e_value=$2;		shift 2 ;;
		--outdir)		outdir=$2;		shift 2 ;;
		--)		break;;
		*)		echo "Unknown option $1 $2" >&2 ; exit 1 ;;
	esac
done

############################################################################################################
[ -z $seq ] && seq="/home/sswang/project/microRNA/sequence/miRNA_seq/miRNA_ATH_miRBase_v20.stem_loop.seq"
[ -z $search_4_dupli_and_single_prog ] && search_4_dupli_and_single_prog="/home/sswang/tools/self_bao_cun/search_4_dupli_and_single/search_4_dupli_and_single.py"
[ -z $blast_result_dir ] && echo "blast_result_dir has to be specified" && exit
[ -z $outdir ] && outdir="miRNA_homolog_result"
[ -d $outdir ] && echo "outdir $outdir has already exists. Exiting ......" && exit 

for i in $blast_result_dir/*blast; do
	echo $i
	basename=`basename $i`
	outname=${basename%%.*blast}
	cmd="
	python $search_4_dupli_and_single_prog \
		--seq_file=$seq \
		--blast_file=$i \
		--outdir=$outdir/$outname \
		--duplicate_evalue=$e_value --singleton_evalue=$e_value \
		--min_coverage=0.5 --min_bit_score=50 --seq_file_format=fasta \
		--coverage_query \
		--force
	"
	${cmd}
done


