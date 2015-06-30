#! /usr/bin/perl
BEGIN{
	my $dir=`dirname $0`;
    chomp($dir);
	push @INC, "$dir/perl_modules";
}

#######################################################################
use 5.010;
use Getopt::Long;
use micro;
use Chang_Yong;
use read_ATH_locus;
use strict;
no autovivification;

my ($psRNAtarget_file, $psRNAtarget_5_file, $UEA_file, $tapir_file, @exclude_type);
my ($exp_file, $outdir, $core_miRNA_swi);
my ($target_result_dir, $target_result_dir, $longest_list_file);
my ($longest_list, $longest_rela_href, $longest_fullname_list_href);
my (@para_type, $para_list_file_dir, $WGD_list_file, $SSD_list_file, $TD_list_file, %para_list);
my (%pair, %para, $para_href, $pair_href);
my (
	$psRNAtarget_result_href, %psRNAtarget_result,
	$UEA_result_href, %UEA_result,
	$tapir_result_href, %tapir_result,
	%predict, %predict_target, $predict
	);
my ($miRNA_psRNAtarget_href, %miRNA_psRNAtarget, $exp_href, %exp, %miRNA_psRNAtarget_5);
my (%locus, $symbol_info_href);
my ($OUT_PAIR_YES, $OUT_PAIR_NO);
my $core_swi=1;
my $force;
my (@compare_type, @duplicate_type, $is_exclude_psRNAtarget_5, $psRNAtarget_5_cutoff);

GetOptions(
	'target_result_dir=s'	    =>	\$target_result_dir,
	'psRNAtarget_file=s'	    =>	\$psRNAtarget_file,
	'exp_file=s'		        =>	\$exp_file,
	'no=s'			            =>	\@exclude_type,
	'outdir=s'          		=>	\$outdir,
	'core_miRNA_swi!'       	=>	\$core_miRNA_swi,
	'force!'	            	=>	\$force,
	'compare_type=s'	        =>	\@compare_type,
	'duplicate_type=s'	        =>	\@duplicate_type,
	'exclude_psRNAtarget_5!'	=>	\$is_exclude_psRNAtarget_5,
	'psRNAtarget_5_cutoff=s'	=>	\$psRNAtarget_5_cutoff,
) || die "illegal param!";

if (not $outdir){
	$outdir = "exp_result.".$$ and `mkdir -p $outdir`;
}
else{
	if (! -d $outdir){`mkdir -p $outdir`;}
	else{
		if ($force){ `rm -rf $outdir`; `mkdir -p $outdir`; }
		else{ die "outdir $outdir has already existed. Exiting ......"; }
	}
}

if (not @compare_type){
	die "compare_type has to be specified! Exiting ......";
}

if (not @duplicate_type){
	die "duplicate_type has to be specified! Exiting ......";
}
else
{
	if ($duplicate_type[0] eq 'all'){@duplicate_type[0..2] = qw(WGD SSD TD);}
}

#################################################################################
#	==================================================================	#
$target_result_dir =	"/home/sswang/project/microRNA/official/target_result" if not $target_result_dir;
$psRNAtarget_file  =	"$target_result_dir/psRNATarget/psRNATargetJob_CDS.txt" if not $psRNAtarget_file;
$psRNAtarget_5_file=	"$target_result_dir/psRNATarget/psRNATargetJob-5.0.txt" if not $psRNAtarget_5_file;
$UEA_file =		"$target_result_dir/UEA_sRNA/UEA_result.txt";
$tapir_file =		"$target_result_dir/tapir/tapir_result_3.5.txt";

$longest_list_file =	"/home/sswang/project/microRNA/sequence/seq_list/TAIR9_longest.list";
$para_list_file_dir=	"/home/sswang/project/microRNA/sequence/seq_list/" if not $para_list_file_dir;
#	==================================================================	#

#open ($OUT_PAIR_YES, '>' ,"YES.list") || die;
#open ($OUT_PAIR_NO,  '>' ,"NO.list")  || die;

{
	@para_type=('WGD','SSD','TD');
	$WGD_list_file = "$para_list_file_dir/WGD_wolfe_pair" if not $WGD_list_file;
	#$WGD_list_file = "$para_list_file_dir/WGD.freeling_2006.list";
	$SSD_list_file = "$para_list_file_dir/SSD_my_new_pair" if not $SSD_list_file;
	$TD_list_file  = "$para_list_file_dir/TD_mayer_pair" if not $TD_list_file;
	%para_list=(
		'WGD' => $WGD_list_file,
		'SSD' => $SSD_list_file,
		'TD'  => $TD_list_file,
	);
	my %duplicate_type_hash_tmp;
	@duplicate_type_hash_tmp{@duplicate_type} = (1) x 10;
	map {print $para_list{$_}."\n"} keys %para_list;
	map {delete $para_list{$_} if not exists $duplicate_type_hash_tmp{$_}} keys %para_list;
	map {print $para_list{$_}."\n"} keys %para_list;
}

#######################################################################################
$symbol_info_href = &read_ATH_locus();

($longest_rela_href, $longest_fullname_list_href) = &obtain_longest($longest_list_file);
	map {$locus{$_}{'longest'}=$longest_rela_href->{$_}} keys %$longest_rela_href;

foreach my $name (keys %para_list){
	my $list_file = $para_list{$name};	# $name stands for 'WGD','SSD','TD'
	($para_href, $pair_href) = &read_para_list($list_file, $name);
	%{$pair{$name}} = %$pair_href;
}

$exp_href = &get_miRNA_target_rela_exp($exp_file);
$predict_target{'exp'}		=	&get_miRNA_target_rela_exp($exp_file);
#foreach (keys %{$predict_target{'exp'}}){print $_."\n";}
$predict_target{'psRNAtarget'}	=	&get_psRNAtarget($psRNAtarget_file, $core_swi);
#$predict_target{'psRNAtarget_5'}=	&get_psRNAtarget($psRNAtarget_5_file, $core_swi, $psRNAtarget_5_cutoff);
$predict_target{'UEA'}		=	&get_UEA($UEA_file, $core_swi);
$predict_target{'tapir'}	=	&get_tapir($tapir_file, $core_swi);

$miRNA_psRNAtarget_href = &get_miRNA_target_rala_psRNAtarget(\%psRNAtarget_result);
	%miRNA_psRNAtarget = %$miRNA_psRNAtarget_href;

#	------------------------------------------------	#
foreach my $symbol (keys %$exp_href){
	my ($gene);
	my $gene = $locus{$symbol}{'longest'};
	if (not $symbol){next;}
	if (not exists $miRNA_psRNAtarget{$symbol}){
		#print "not exists $symbol\n";
		next;
	}
	for (keys %{$exp_href->{$symbol}}){
		$_ = &del_taxaname_miRNA($_);
		$_ = &get_core_miRNA($_) if $core_miRNA_swi;
		#print $symbol."\t".$_."\n" if not exists $miRNA_psRNAtarget{$symbol}{$_};
	}
}

&get_Two_Outof_Three(\%predict_target, $outdir);

{
	open (my $OUT_both, '>', "$outdir/both.list");
	open (my $OUT_either, '>', "$outdir/either.list");
	open (my $OUT_neither, '>', "$outdir/neither.list");
    my %num_of_targets_in_a_pair;
    my $num_of_targets_in_a_pair_href=\%num_of_targets_in_a_pair;
	foreach my $type (keys %pair){
		open (my $OUT_both_type, '>', "$outdir/both_${type}.list");
		open (my $OUT_either_type, '>', "$outdir/either_${type}.list");
		open (my $OUT_neither_type, '>', "$outdir/neither_${type}.list");
		foreach my $pair (keys %{$pair{$type}}){
			my @pair = split /\-/, $pair;
			#&generate_exp($pair, $exp_href);
			#&generate_predict($pair, ['2_outof_3'], [$OUT_both, $OUT_either, $OUT_neither]);
			$num_of_targets_in_a_pair_href = &generate_predict($pair, \@compare_type, [$OUT_both, $OUT_either, $OUT_neither],\%num_of_targets_in_a_pair);
			&generate_predict($pair, \@compare_type, [$OUT_both_type, $OUT_either_type, $OUT_neither_type],\%num_of_targets_in_a_pair);
            #print $num_of_targets_in_a_pair_href->{$type}{'1'}."\n";
		}
		close $OUT_both_type; close $OUT_either_type; close $OUT_neither_type;
		print STDOUT $type . "\t" . scalar(keys %{$pair{$type}}) . "\n";
	}
        for my $type (keys %$num_of_targets_in_a_pair_href){
            print "##"."\t".$type."\n";
            for my $num (keys %{$num_of_targets_in_a_pair_href->{$type}}){
                print STDOUT $type."\t".$num."\t".$num_of_targets_in_a_pair_href->{$type}{$num}."\n";
            }
        }
	close $OUT_both; close $OUT_either; close $OUT_neither;
}

#######################################################################################
sub generate_predict{
	my ($pair, $HashName_aref, $OUT_aref, $num_of_targets_in_a_pair_href) = @_;
	my ($OUT_both, $OUT_either, $OUT_neither) = @$OUT_aref;
	my (%divergent_miRNA_hash, %all_miRNAs);
	my @pair = split /-/, $pair;
	my $pair_k=0; # if 有一个paralog有实验证据就+1
	my $pair_divergence=0; # if 有一个divergent miRNA binding sites pattern 就+1
	my $pair_no_both=0; # HashName必须相同，这里排除了psRNAtarget_5的干扰，也就是说如果没在psRNAtarget_5中出现但是在比如'2_outof_3'中出现了就可能是both
	my @HashName = @$HashName_aref;
    my $num_of_targets=0;
	#my @HashName = qw(psRNAtarget UEA tapir);
	for my $HashName (@HashName){
		for my $order (0..$#pair){
            my $symbol  = $pair[$order];
			my $paralog = $pair[1-$order];
			if (exists $predict_target{$HashName}{$symbol}){
                $num_of_targets++; # if %{$predict_target{$HashName}{$symbol}};
				foreach my $miRNA (keys %{$predict_target{$HashName}->{$symbol}}){
					$all_miRNAs{$miRNA}=1;
					$pair_k++;
					my $pair_divergence_tmp_k=0;
					my @not_existing_HashName = ($HashName); # 这里要注意，需要灵敏调整随机应变
					if ($is_exclude_psRNAtarget_5){
						push @not_existing_HashName, 'psRNAtarget_5' if $HashName ne 'exp';
					}
					for my $name (@not_existing_HashName){
						if (not exists $predict_target{$name}->{$paralog}{$miRNA}){
							$divergent_miRNA_hash{$miRNA}=1;
							$pair_divergence_tmp_k++;
							$pair_divergence++;
							$pair_no_both++ if $name eq $HashName;
						}
					}
					#$pair_divergence++ if $pair_divergence_tmp_k != scalar(@not_existing_HashName);
					#print $pair_divergence."\t".scalar(@not_existing_HashName)."\n";
				}
			}
		}
	}

	if ($pair_k != 0){
		if ($pair_divergence>=1){
			print $OUT_either $pair . "\t";
			print $OUT_either join("\t", keys %divergent_miRNA_hash) . "\n";
            $num_of_targets_in_a_pair_href->{"hashname"}{$num_of_targets}++;
            print $pair."\n" if $num_of_targets==2;
		}
		if ($pair_divergence == 0){
			print $OUT_both $pair . "\t";
			print $OUT_both join("\t", keys %all_miRNAs) . "\n";
		}
	}
	else{
		print $OUT_neither $pair."\n";
	}

    return ($num_of_targets_in_a_pair_href);
}


sub get_Two_Outof_Three{
	my (%two_outof_three, %filter_target);
	my ($predict_target_href, $outdir) = @_;
	open (my $OUT_2_outof_3, '>', "$outdir/target.2_outof_3.list");
	for my $type (keys %$predict_target_href){
		open (my $OUT, '>', "$outdir/target.$type.list");
		next if $type eq 'psRNAtarget_5'; # exclude psRNAtarget_5 when counting
		print $type."\n";
		next if $type eq 'exp';
		foreach my $symbol (sort keys %{$predict_target_href->{$type}}){
			print $OUT join("\t", $symbol, sort keys %{$predict_target_href->{$type}{$symbol}});
			print $OUT "\n";
			for my $miRNA (keys %{$predict_target_href->{$type}{$symbol}}){
				$filter_target{$symbol}{$miRNA} ++;
			}
		}
		close $OUT;
	}
	for my $symbol (sort keys %filter_target){
                my @miRNAs;
		for my $miRNA (sort keys %{$filter_target{$symbol}}){
			if ($filter_target{$symbol}{$miRNA} >= 2){
				$predict_target{'2_outof_3'}{$symbol}{$miRNA} = 1;
				@miRNAs = keys %{$predict_target{'2_outof_3'}{$symbol}};
			}
		}
		print $OUT_2_outof_3 join ("\t", $symbol, @miRNAs), "\n" if @miRNAs;
		#my @miRNA = keys %{$predict_target{'2_outof_3'}{$symbol}};
		#print $OUT join ("\t", $symbol, @miRNA), "\n" if @miRNA;
	}
	close $OUT_2_outof_3;
	my $num = scalar (keys %{$predict_target{'2_outof_3'}}); print "2_outof_3 num is ".$num."\n";
}

sub generate_exp{
	my $pair_k; # if 有一个paralog有实验证据就+1
	my ($pair, $exp_href) = @_;
	my @pair = split /-/, $pair;
	for my $order (0..$#pair){
		my $symbol  = $pair[$order];
		my $paralog = $pair[1-$order];
		if (exists $exp_href->{$symbol}){
			$pair_k++;
			#print "$pair\n" if $pair_k == 1;
			foreach (keys %{$exp_href->{$symbol}}){
				my ($print_out);
				my $miRNA = $_;
				$print_out = $miRNA."\t";
				map {$print_out.=$_."\t".$symbol_info_href->{$_}{'description'}."\t###\t"} ($symbol, $paralog);
				$print_out .= "\t";
				if (exists $miRNA_psRNAtarget_href->{$paralog}{$miRNA} or exists $exp_href->{$paralog}{$miRNA}){
					$print_out .= "YES"
				}
				else{
					$print_out .= "NO";
				}
				$print_out =~ /\t([^\t]+)$/;
=cut
				given ($1){
					when ($1 eq 'YES')	{print $OUT_PAIR_YES $pair."\n"}
					when ($1 eq 'NO')	{print $OUT_PAIR_NO  $pair."\n"}
				}
=cut
			}
		}
	}
}


sub get_miRNA_target_rela_exp{
	my %exp;
	my ($exp_file) = @_;
	open (my $IN, '<', "$exp_file") || die "exp_file cannot be opened";
	while(<$IN>){
		chomp;
		my ($symbol);
		#my ($miRNA, $symbol) = split /\t/;
		my ($symbol, @miRNA) = split /\t/;
		for my $miRNA (@miRNA){
			$miRNA = &del_taxaname_miRNA($miRNA);	# ath-miR1234 -> miR1234
			$miRNA = &get_core_miRNA($miRNA) if $core_miRNA_swi;
			$exp{$symbol}{$miRNA} = 1;
			$predict_target{'exp'}{$symbol}{$miRNA}=1;
		}
		#$exp{$symbol}{$miRNA} = 1;
	}
	return (\%exp);
}

sub get_miRNA_target_rala_psRNAtarget{
	my (%miRNA_psRNAtarget_5);
	my ($ref_method)=@_;
	foreach (keys %$ref_method){
                my $miRNA=$ref_method->{$_}{'miRNA'};
		$_ = &get_core_name($_);
		$miRNA = &get_core_miRNA($miRNA);
		#print $miRNA."\n";
		$miRNA_psRNAtarget_5{$_}{$miRNA}++;
        }
	return (\%miRNA_psRNAtarget_5);
}


# ======================= Get Original Prediction Results =========================== #
# =================================================================================== #
sub get_psRNAtarget
{
my ($psRNAtarget_file, $core_swi, $score_cutoff) = @_;
$score_cutoff=10 if not $score_cutoff;
my ($count_line,%psRNAtarget,%psRNAtarget_result);
print $psRNAtarget_file."\n";
open(my $IN , '<' , "$psRNAtarget_file") || die "psRNAtarget_file $psRNAtarget_file cannot be opened";
while(<$IN>){
	next if /^#/;
    next if /^miRNA_Acc\.\t/;
    chomp();
	#$_ =~ s/\r?\n//;
	my ($miRNA,$seq_name,$score) = split /\t/;
	next if $score >= $score_cutoff;
	($seq_name) = &uc($seq_name);
	$miRNA = &del_taxaname_miRNA($miRNA);	# ath-miR1234 -> miR1234
	$miRNA = &get_core_miRNA($miRNA) if $core_miRNA_swi;
	$seq_name = &get_core_name($seq_name) if $core_swi;
	$psRNAtarget_result{$seq_name}{$miRNA} = 1;
	#$psRNAtarget_result{$seq_name}{'score'} = $score;
}
close $IN;
return (\%psRNAtarget_result);
}

sub get_UEA
{
my ($UEA_file, $core_swi) = @_;
open(IN,'<',"$UEA_file") || die "UEA input file $UEA_file cannot be opened!";
while(<IN>){
	chomp;
	my ($miRNA,$seq_name,$posi) = split /\,/;
    next if $posi eq "No target found";
	# $posi is the posi where target site is predicted to be
	my ($seq_name) = &uc($seq_name);
	$miRNA = &del_taxaname_miRNA($miRNA);   # ath-miR1234 -> miR1234
	$miRNA = &get_core_miRNA($miRNA) if $core_miRNA_swi;
	$seq_name = &get_core_name($seq_name) if $core_swi;
	$UEA_result{$seq_name}{$miRNA} = 1;
	#$UEA_result{$seq_name}{'posi'} =$posi;
}
close IN;
return (\%UEA_result);
}

sub get_tapir
{
my (%tapiur_result);
my ($tapir_file, $core_name) = @_;
open (my $IN_TAPIR,'<',$tapir_file) || die "tapir input file $tapir_file cannot be opened!";
while(<$IN_TAPIR>){
	chomp;
	my ($miRNA, $seq_name, $score);
	next if not /^miRNA.+(miR\d+)/;
	$miRNA=$1;
	$miRNA = &get_core_miRNA($miRNA) if $core_miRNA_swi;
	while(<$IN_TAPIR>){
		chomp;
        $seq_name=$1 if $_ =~ /^target\s+(.+)\S?/;
		($seq_name) = &uc($seq_name);
        $score=$1 and last if $_ =~ /score\s+(.+)/;
    }
	$miRNA = &del_taxaname_miRNA($miRNA);   # ath-miR1234 -> miR1234
	$seq_name = &get_core_name($seq_name) if $core_swi;
	$tapir_result{$seq_name}{$miRNA}=1;
	#$tapir_result{$seq_name}{'score'}=$score;
}
close $IN_TAPIR;
return (\%tapir_result);
}


