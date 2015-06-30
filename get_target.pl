#! /usr/bin/perl
use 5.010;
use Getopt::Long;

my $longest_swi='';
my ($target_result_dir,$psRNAtarget_file,$tapir_file,$UEA_file,$psRNAtarget_5_file,@target_result);
my ($psRNAtarget_result);
my (%para_list,$para_list_file_dir,$WGD_list_file,$SSD_list_file,$TD_list_file,$WGD,$SSD,$TD,@para_type);
my ($longest_list_file,$ref_longest_rela,$ref_longest_full_name);
my (%miRNA);

##################################################################
GetOptions(
	'longest'		=>	\$longest_swi,
	'target_result_dir'	=>	\$target_result_dir,
	'para_list_file_dir'	=>	\$para_list_file_dir,
) || die "illegal options!";

##################################################################

my $target_list_output = 'all_target.list';
open (my $TARGET_OUT, '>', $target_list_output);

$longest_list_file="../sequence/seq_list/TAIR9_longest.list";

{
$target_result_dir = "link/official/target_result";
$psRNAtarget_file = "$target_result_dir/psRNATarget/psRNATargetJob_CDS.txt";
$psRNAtarget_5_file = "$target_result_dir/psRNATarget/psRNATargetJob-5.0.txt";
$tapir_file = "$target_result_dir/tapir/tapir_result_3.5.txt";
$UEA_file = "$target_result_dir/UEA_sRNA/UEA_result.txt";
}

{
	@para_type=('WGD','SSD','TD');
	$para_list_file_dir="../sequence/seq_list/";
	$WGD_list_file = "$para_list_file_dir/WGD_wolfe_pair";
	$SSD_list_file = "$para_list_file_dir/SSD_mayer_pair";
	$TD_list_file  = "$para_list_file_dir/TD_mayer_pair";
	%para_list=(
		'WGD' => $WGD_list_file,
		'SSD' => $SSD_list_file,
		'TD'  => $TD_list_file,
	);
}

##########################################################
foreach my $name(keys %para_list){
	my $list_file = $para_list{$name};	# $name stands for 'WGD','SSD','TD'
	${'ref_para_'.$name} = &read_para_list($list_file, $name);
}

{
	*psRNAtarget_result	= &get_psRNAtarget($psRNAtarget_file);
	*psRNAtarget_5_result = &get_psRNAtarget($psRNAtarget_5_file);
	*UEA_result = &get_UEA($UEA_file);
	*tapir_result = &get_tapir($tapir_file);
}


($ref_longest_rela,$ref_longest_full_name) = &generate_longest($longest_list_file);	#foreach (keys %$ref_longest_rela) {print $_."\n";}

#*miRNA = &get_miRNA_target_rela(*psRNAtarget_result,*UEA_result,*tapir_result);
(*miRNA,*miRNA_UEA,*miRNA_psRNAtarget,*miRNA_tapir) = &get_miRNA_target_rela(*psRNAtarget_result,*UEA_result,*tapir_result);
*miRNA_psRNAtarget_5 = &get_miRNA_target_rala_psRNAtarget_5(*psRNAtarget_5_result);

*target = &get_target(\%miRNA, 2, 100, $longest_swi);
*target_psRNAtarget_5 = &get_target(\%miRNA_psRNAtarget_5, 1, 1, $longest_swi);

#foreach(keys target_psRNAtarget_5){print $_."\n"}


foreach my $seq_name(keys %{$target{'2-100'}}){
	next if not exists $para_WGD{$seq_name};
	my $para_name = $ref_para_WGD->{$seq_name};
	next if exists $target{'2-100'}{$para_name};
	#print $seq_name."\t".$para_name."\n";
	#print $seq_name."\n" if exists ${'para_WGD'}{$seq_name} and not exists ${'para_WGD'}{$para};
}

&PRINT_target_list(\%{$target{'2-100'}}, $TARGET_OUT);

foreach (@para_type){
	generate_result($_);
}

&generate_final_result(\%final_num, );

##############################################################
##############################################################
sub generate_result
{
my ($type_duplication) = @_;
my $pair_hashname = 'pair_'.$type_duplication;
foreach (keys %{$pair_hashname}){
	my (%count, @symbol, $num_keys_count, $count_order);
	# count_order is the order (1 or 2) or the symbok targeted by microRNA(s) in the array @symbol
	@symbol =  split /\-/;
	for (0..$#symbol){
		$count{$symbol[$_]}=1 and $count_order=$_ if exists $target{'2-100'}{$symbol[$_]};
	}
	
	$num_keys_count = scalar keys %count;
	if ($num_keys_count == 2){
		$final_num{$type_duplication}{both_with_divergent}++;
		#print "@symbol\n";
		my (@miRNA_tmp);
		for my $symbol (@symbol){
			push @miRNA_tmp, \%{$miRNA{$symbol}};
		}
		if (%{$miRNA_tmp[0]} ~~ %{$miRNA_tmp[1]}){
			#print "@symbol\n";
			$final_num{$type_duplication}{both_no_divergent}++;
		}
		else{
			$final_num{$type_duplication}{both_divergent}++;
		}
	}

	elsif ($num_keys_count == 1){
		print $symbol[$count_order]."\t".$symbol[1-$count_order]."\t".$type_duplication."\n";
		$final_num{$type_duplication}{either}++;
		$final_num{$type_duplication}{either_strict}++ if not exists $target_psRNAtarget_5{'1-1'}{$symbol[$count_order]};
	}
	
	else{
		$final_num{$type_duplication}{neither}++;
	}
}
}

sub generate_final_result
{
my ($final_num_hashref)=@_;
my %final_num = %$final_num_hashref;
foreach my $type_duplication (keys %final_num){
	print $type_duplication."\n";
	foreach (keys %{$final_num{$type_duplication}}){
		print "$_\t".$final_num{$type_duplication}{$_}."\n";
	}
	print "\n";
}
}

sub PRINT_target_list{
	my ($hashref, $outfile_FH)=@_;
	select $outfile_FH;
	for my $symbol (keys %$hashref){
		print $symbol."\t";
		print join "\t", keys %{$miRNA{$symbol}};
		print "\n";
	}
	select STDOUT;
}

sub get_psRNAtarget
{
my ($psRNAtarget_file) = @_;
my ($count_line,%psRNAtarget,%psRNAtarget_result);
open(my $IN , '<' , $psRNAtarget_file);
while(<$IN>)
{
	next if ++$count_line <= 2;
	$_ =~ s/\r?\n//;
	my ($miRNA,$seq_name,$score) = split /\t/;
	&uc($seq_name);
	$miRNA = &del_taxaname_miRNA($miRNA);	# ath-miR1234 -> miR1234
	$psRNAtarget_result{$seq_name}{'miRNA'} = $miRNA;
	$psRNAtarget_result{$seq_name}{'score'} = $score;
}
close $IN;
return (\%psRNAtarget_result);
}


sub get_UEA
{
my ($UEA_file) = @_;
open(IN,'<',"$UEA_file");
while(<IN>){
	chomp;
	my ($miRNA,$seq_name,$posi) = split /\,/;
	# $posi is the posi where target site is predicted to be
	&uc($seq_name);
	$miRNA = &del_taxaname_miRNA($miRNA);   # ath-miR1234 -> miR1234
	$UEA_result{$seq_name}{'miRNA'}=$miRNA;
	$UEA_result{$seq_name}{'posi'} =$posi;
}
close IN;
return (\%UEA_result);
}


sub get_tapir
{
my ($tapir_file) = @_;
open (my $IN_TAPIR,'<',$tapir_file);
while(<$IN_TAPIR>){
	chomp;
	my ($miRNA,$seq_name);
	next if not /^miRNA.+(miR\d+)/;
	$miRNA=$1;
	while(<$IN_TAPIR>)
       	{
		chomp;
              	$seq_name=$1 if $_ =~ /^target\s+(.+)\S?/;
              	$score=$1 and last if $_ =~ /score\s+(.+)/;
       	}
	$miRNA = &del_taxaname_miRNA($miRNA);   # ath-miR1234 -> miR1234
	$tapir_result{$seq_name}{'miRNA'}=$miRNA;
	$tapir_result{$seq_name}{'score'}=$score;
}
close $IN_TAPIR;
return (\%tapir_result);
}


sub read_para_list{
	my ($list_file,$list_name) = @_;
	#my %{'para_'.$list_name};
	open (my $IN, '<', $list_file);
	while(<$IN>){
		chomp;
		$_=uc;
		my @line=split;
		next if scalar (@line) >= 3;
		foreach (0..1){
			next if not defined $line[$_];
			${'para_'.$list_name}{$line[$_]}=$line[1-$_];
		}
		${'pair_'.$list_name}{join "-",sort @line[0,1]} = 1;
		#${$list_name}{join '-',sort split}=1;
	}
	return (\%{'para_'.$list_name}, \%{'pair_'.$list_name});
}


sub generate_longest{
	my ($list_file) = @_;
	my (%longest_rela,%longest_full_name);
	open(my $IN,'<',$list_file);
	while(<$IN>){
		my ($core_name);
		chomp;
		($core_name) = $_ =~ /([^\.]+)\.\d+/;
		$longest_rela{$core_name}=$_;
		$longest_full_name{$_}=1;
	}
	return(\%longest_rela,\%longest_full_name);
}


sub get_miRNA_target_rela
{
my (@lala);
my @target_result = @_;
foreach my $ref_method(sort @target_result){
	# $ref_method stands for the REF of prediction method used, i.e. tapir
	foreach (keys %$ref_method){
		my ($core_name, $miRNA);
		$core_name = get_core_name($_);
		$miRNA=$ref_method->{$_}{'miRNA'};
		$miRNA{$_}{$miRNA}++;
		$miRNA{$core_name}{$miRNA}++;
		${'miRNA_'.$ref_method}{$_}{$miRNA}++;
	}
}
return(\%miRNA, do{
	foreach(sort @target_result) {push @lala,\%{'miRNA_'.$_};}
	@lala;
	}
);
}


sub get_target
{
my (%target);
my ($REF_miRNA,$lower_limit,$upper_limit) = @_;
foreach my $seq_name(keys %$REF_miRNA){
	#next if not exists $ref_longest_full_name->{$seq_name};
	my $core_name = &get_core_name($seq_name);
	foreach my $miRNA(keys %{${$REF_miRNA}{$seq_name}}){
		my $miRNA_num=${$REF_miRNA}{$seq_name}{$miRNA};
		#$target{$lower_limit.'_'.$upper_limit}{$core_name}=1 if ($miRNA_num >= $lower_limit and $miRNA_num <= $upper_limit);
		$target{$lower_limit.'-'.$upper_limit}{$core_name}=1 if ($miRNA_num >= $lower_limit and $miRNA_num <= $upper_limit);
	}
}
return (\%target);
}

##################################################################
sub uc{
	foreach (@_) {$_=uc($_)};
}

sub del_taxaname_miRNA   # ath-miR1234 -> miR1234
{
	my ($miRNA_name) = @_;
	$miRNA_name =~ s/^[^-]+\-(?=miR\d+)//;
	return ($miRNA_name);
}

sub get_core_name{
	my ($full_name) = @_;
	my ($core_name) = $full_name =~ /([^\.]+)\.?/;
	return ($core_name);
}

sub get_core_miRNA{
	$_[0] =~ s/\D+$//;
	return ($_[0]);
}

#-------------------------------------------------------------#
#----------      not very important func		------#
sub get_miRNA_target_rala_psRNAtarget_5{
	my ($ref_method)=@_;
	foreach (keys %$ref_method){
                my $miRNA=$ref_method->{$_}{'miRNA'};
                $miRNA_psRNAtarget_5{$_}{$miRNA}++;
        }
	return (\%miRNA_psRNAtarget_5);
}


##################################################################
if (1<0){
foreach (keys %$ref_psRNAtarget_result)	{print $_."\n"}
foreach (keys %$ref_UEA_result)		{print $_."\n"}
foreach (keys %$ref_tapir_result)	{print $_."\n"}
}


####################################################################
sub generate_result_old
{
my ($type_duplication) = @_;
foreach my $seq_name(keys %{$target{'2-100'}}){
	my $core_name=&get_core_name($seq_name);
	#print $core_name."\n";
	#next;
	foreach my $list_name(@para_type){
		next if $list_name ne $type_duplication;
		my $para_type='para_'.$list_name;
		next if not exists ${$para_type}{$core_name};
		my $para_seq_name=${'ref_'.$para_type}->{$core_name};
		#print $para_seq_name."\n";
		next if not exists $target{'2-100'}{$para_seq_name};

		$para_core_name=&get_core_name($para_seq_name);
		my $a=join "\t",sort ($para_core_name,$seq_name);
		#print $para_type."\t".$a."\n";
		print $a."\n";
		next;
		if (not exists $target{'1-100'}{$para_core_name}){
			print $para_type."\t".$seq_name."\n";
		}
	}
}
}

