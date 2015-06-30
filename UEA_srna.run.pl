#! /usr/bin/perl

use strict;
use Getopt::Long;
use File::Basename;

my ($srna_file, $gene_file, $UEA_PATH, $out_dir, $output);
my ($new_gene_file);
#############################################

GetOptions(
	'srna_file=s'	=>	\$srna_file,
	'gene_file=s'	=>	\$gene_file,
	'UEA_PATH=s'	=>	\$UEA_PATH,
	'out_dir=s'	=>	\$out_dir,
	'output=s'	=>	\$output,
) || die "illegal param";

$UEA_PATH="/home/sswang/software/RNA/srna-2/program" if not $UEA_PATH;
	die "UEA_PATH does not exist" if ! -d $UEA_PATH;
die "mandatory files not defined" if not $gene_file or not $srna_file;

$out_dir = "UEA_outdir" if not $out_dir;
	`rm -rf "$out_dir"` if -e $out_dir;
	mkdir "$out_dir" if not -e $out_dir;

system "cp \"$gene_file\" \"$UEA_PATH/data/transcriptomes\"";
$new_gene_file = basename("$gene_file");
$output = "$new_gene_file".'.UEA_output';

################################################

system "\"$UEA_PATH\"/srna-tools.pl --tool target --out \"$out_dir/$output\" --srna_file \"$gene_file\" --transcriptome \"$new_gene_file\"";

`sed -i '/^[0-9].*No target found/d' \"$out_dir/$output/results.csv\"`;



