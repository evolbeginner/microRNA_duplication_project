#! /bin/bash

#####################################################################
blastn=/home/sswang/bin/blastn
blastdb_dir="/mnt/storage2/alex/liu_linc/data/blast_dbs";
microRNA_seq_file="/home/sswang/project/microRNA/sequence/miRNA_seq/miRNA_ATH_miRBase_v20.stem_loop.seq"
cpu=2
e_value=1e-10
outdir=`pwd`/blast_result

#####################################################################
function get_Params(){
	getopt_cmd="getopt -o '' --long help,blastdb_dir:,e_value:,cpu:,CPU:,outdir:,microRNA_seq_file:,force, -n 'bmtagger' -- "$@""
	PARAM=$($getopt_cmd)
	eval set -- "$PARAM"
	while true ; do
		case "$1" in
			#-h|--help)	show_help;		shift ;;
			--blastdb_dir)	blastdb_dir=$2;		shift 2 ;;
			--e_value)	e_value=$2;		shift 2;;
			--cpu|--CPU)	cpu=$2;			shift 2;;
			--outdir)	outdir=$2;		shift 2;;
			--microRNA_seq_file) microRNA_seq_file=$2; shift 2;;
			--force)	force=1;		shift;;
			--)		break;;
			*)		echo "Unknown option $1 $2" >&2 ; exit 1 ;;
		esac
	done

	if [ -d $outdir ]; then
		if [ ! -z $force ]; then
			rm -rf $outdir
			mkdir $outdir
		fi
	else
		mkdir $outdir
	fi
}

function run_Blast(){
	declare -A db_name_set
	for i in `ls $blastdb_dir`; do
		#for nin,nhr,nsq
		db_name=${i%%\.*}
		db_name_set[$db_name]=1
		microRNA_seq_file_basename=`basename $microRNA_seq_file`
		$blastn -query $microRNA_seq_file -db $blastdb_dir/$db_name -outfmt 6 -evalue $e_value -num_threads $cpu -out $outdir/${db_name}.blast
	done
}

#####################################################################

get_Params $@

run_Blast


