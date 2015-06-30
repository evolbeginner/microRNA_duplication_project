#! /usr/bin/perl

use Getopt::Long;
use Statistics::Basic qw(:all);
use Statistics::PointEstimation;
use List::Util qw(sum first);
use Statistics::TTest;
use Statistics::Test::WilcoxonRankSum;
use File::Basename;

use 5.010;


##################################################################################
my ($seq_list1, $seq_list2, $seq_list_file_separator, @ka_file, $ka_cutoff, $ks_cutoff, $ka_swi, $ks_swi, $kaks_swi,
	$rbs_swi, $range_swi, $range_upper, $range_step, $force_swi,
	$seq_similarity, $seq_similarity_cutoff, $seq_similarity_file,
    $is_detail, $outdir, $is_force, $is_auto_outfile, $outfile1, $outfile2);
my ($rbs_file);
my (%seq_similarity_pass);

$ka_swi = True;
$is_auto_outfile = True;


GetOptions(
    'in1=s' 			        =>	\$seq_list1,
    'in2=s' 			        =>	\$seq_list2,
    'seq_list_file_separator=s' =>	\$seq_list_file_separator,
    'ka_file=s'			        =>	\@ka_file_arr,
    'ka_cutoff=s'			    =>	\$ka_cutoff,
    'ks_cutoff=s'			    =>	\$ks_cutoff,
    'ka_swi!'			        =>	\$ka_swi,
    'ks_swi!'			        =>	\$ks_swi,
    'kaks_swi!'			        =>	\$kaks_swi,
    'rbs_swi!'                  =>  \$rbs_swi,
    'range=s'			        =>	\$range_swi,
    'range_upper=s'			    =>	\$range_upper,
    'range_step=s'			    =>	\$range_step,
    'range_freq=s'			    =>	\$range_freq,
    'force!'			        =>	\$force_swi,
    'seq_similarity!'		    =>	\$seq_similarity,
    'seq_similarity_cutoff=s'	=>	\$seq_similarity_cutoff,
    'seq_similarity_file=s'		=>	\$seq_similarity_file,
    'is_split_tab_first!'		=>	\$is_split_tab_first,
    'detail|details!'           =>  \$is_detail,
    'outdir=s'                  =>  \$outdir,
    'force'                     =>  \$is_force,
    'auto_outfile'              =>  \$is_auto_outfile,
) || die "illegal param!\n";


&check_param();
*seq_similarity_pass = &generate_seq_similarity_pass($seq_similarity_file, $seq_similarity_cutoff) if $seq_similarity_file;


if (not defined $seq_list1){
	$seq_list1 = "result_inter_Guan_SSD_physical/connect_para_80-more.list";
	#"/home/sswang/project/dimer/seq/seq_list/WGD_PD.list";
}

if (not defined $seq_list2){
	$seq_list2 = "result_inter_Guan_SSD_physical/connect_non_para_80-more.list";
	#"/home/sswang/project/dimer/seq/seq_list/WGD_non_PD.list";
}


@seq_list_file = ($seq_list1,$seq_list2);
$rbs_swi='ON' if defined $rbs_swi;
$rbs_file="/home/sswang/project/dimer/seq/seq_list/rbs.list";

if (not defined @ka_file_arr){
	#$ka_file = "/home/sswang/project/dimer/result/kaks/yeast_WGD.kaks";
	$ka_file_arr = ("/home/sswang/project/dimer/result/kaks/yeast_Guan_SSD.kaks");
}


*rbs = &extract_rbs($rbs_file);
*ka = &read_ka(\@ka_file_arr);
*ka_hash = &read_seq_list(\@seq_list_file, $seq_list_file_separator);


foreach my $key1(keys %ka_hash){
	my $mean = mean(keys %{$ka_hash{$key1}});
	&range_ka ( $key1 , \%{$ka_hash{$key1}} , $range_upper , $range_step , $range_freq );
	my @data = grep {$_ if $_ =~ /^-?[.0-9]+$/} keys %{$ka_hash{$key1}};	
	#foreach(keys %{$ka_hash{$key1}}) {print $_."\n" }
	#push @stat_test, [keys %{$ka_hash{$key1}}];
	push @stat_test, [@data];
	print $key1."\t".$mean."\n";
}


if (defined $range_swi){
open (my $RANGE_OUT, '>', "$range_swi") || die "range_swi $range_swi cannot be opened!\n";
select $RANGE_OUT;
foreach my $key1(sort keys %range){
	#foreach(keys %{$range{$key1}}){print $_."\n"};
	my $num_1 = $range{$key1}{1};
	my $num_2 = $range{$key1}{2};
	my $mean = &mean(split /\-/,$key1);
	print $key1."\t".$mean."\t".$num_1/($num_1+$num_2)."\t".$num_1."\t".$num_2."\n";	
}
close $RANGE_OUT;
select STDOUT;
}


if($is_detail){
    open(my $OUT_fh, '>', $outfile1) || die "";
    print $OUT_fh join("\n", @{$stat_test[0]});
    close $OUT_fh;

    open(my $OUT_fh, '>', $outfile2) || die "";
    print $OUT_fh join("\n", @{$stat_test[1]});
    close $OUT_fh;
}


$ttest = new Statistics::TTest;
$ttest->load_data(@stat_test);
my $prob = $ttest->{t_prob};
print $prob."\n";
#$ttest->output_t_test();

my $wilcox_test = Statistics::Test::WilcoxonRankSum->new();
$wilcox_test->load_data(@stat_test);
my $prob_wilcox = $wilcox_test->probability();
print $prob_wilcox."\n";


##########################################################################
sub generate_seq_similarity_pass{
	my (%seq_similarity_pass);
	my ($seq_similarity_file, $seq_similarity_cutoff) = @_;
	die "seq_similarity_cutoff $seq_similarity_cutoff has not be specified!\n" if not $seq_similarity_cutoff;
	open ($IN, '<', $seq_similarity_file) || die "seq_similarity_file $seq_similarity_file cannot be opened!\n";
	while(<$IN>){
		chomp;
		my ($gene1, $gene2, $seq_similarity) = split;
		if ($seq_similarity >= $seq_similarity_cutoff){
			my $pair = join ("|", sort ($gene1, $gene2));
			$seq_similarity_pass{$pair} = 1;
		}
	}
	close $IN;
	return (\%seq_similarity_pass);
}


sub read_seq_list{
	my ($seq_list_file_aref, $seq_list_file_separator) = @_;
	$seq_list_file_separator = "\\s+" if not $seq_list_file_separator;
	for(0..$#seq_list_file){
		my $num=$_+1;
		open(my $IN,'<',$seq_list_file[$_]) || die "$seq_list_file $seq_list_file[$_] cannot be opened!\n";
		while(<$IN>){
			chomp;
			my @line = split/\t/;
			my $line = $is_split_tab_first ? $line[0] : $_; # ∫‹÷ÿ“™
			my @seq_pair = (split /$seq_list_file_separator/, $line)[0,1];
			my $seq_pair = join ("|", sort @seq_pair);
			my $ka=$ka{$seq_pair};
			$ka_hash{$num}{$ka}=1;
		}
		close $IN;
	}
	return(\%ka_hash);
}


sub read_ka{
	my ($ka_file_arr_aref) = @_;
	my %ka;
	foreach my $ka_file (@$ka_file_arr_aref){
		open(my $IN,'<',$ka_file) || die "ka_file $ka_file cannot be opened!\n";
		while(<$IN>){
			chomp;
			my @a=split;
			if (not $seq_similarity){
				next if $a[4] !~ /^-?[.0-9]+$/;
			}
			if (defined $rbs_swi and $rbs_swi eq 'ON'){
				next if exists $rbs{$a[0]} or exists $rbs{$a[1]};
			}
			if (defined $ka_cutoff){next if $a[2] > $ka_cutoff};
			if (defined $ks_cutoff){next if $a[3] > $ks_cutoff};
			if (defined $seq_similarity_cutoff and defined $seq_similarity_file){next if not exists $seq_similarity_pass{join "|", sort @a[0,1]}};
			my ($ka, $ks, $ka_ks) = @a[2,3,4];
			$order = 2;
			$order = 2 if $ka_swi;
			$order = 3 if $ks_swi;
			$order = 4 if $kaks_swi;
			$ka{join ("|", sort @a[0,1])} = $a[$order];
		}
	}
	return (\%ka);
}


sub range_ka{
	my ($key, $ref, $range_upper, $range_step, $range_freq) = @_;
	$range_step = 0.2	if not $range_step;
	$range_upper = 1	if not $range_upper;
	$range_freq = 200	if not $range_freq;
	for my $key2(keys %$ref){
		#my @range = map {$_/2} 0..8;
		my @range = qw( 0-2 1-2 2-2.5 2.5-3 3-3.5 3.5-3.8 3.8-4
				3.5-4 4-4.2 4.2-4.4 4.4-4.6 4.6-4.8 4.8-5 5-5.3 5-5.5 5.5-6 );
		my @range = map {
					my $a=0+$_*$range_upper/$range_freq;
					my $b1=$a-$range_step;
					my $b2=$a+$range_step;
					$b1.'-'.$b2;
				} 1..$range_freq;

		for my $num(@range){
			my ($num1,$num2)=split /\-/,$num;
			if ($key2>$num1 and $key2<=$num2){
				$range{$num1.'-'.$num2}{$key}++;
				#print $num1.'-'.$num2."\t".$range{$num1.'-'.$num2}{$key1}."\n";
			}
		}
	}
}


sub mean{
	my @a = map {$_=~/^-?[.0-9]+$/ and $_} @_;
	return sum(@a)/@a;
}


sub extract_rbs
{
	my ($rbs_file)=@_;
	my %rbs;
	open(IN,'<',$rbs_file);
	while(<IN>){
		chomp;
		do {      $_=uc($_);	$rbs{$_} = 1    } foreach (split);
	}
	close IN;
	return(\%rbs);
}


sub check_param{
	my $ji_shu;
	for my $i(qw(seq_list1 seq_list2 ka_file)){
		print $i."\n" and $ji_shu++ if not defined \${$i};
	}
	print $ji_shu."\n";
	return if defined $force_swi;
	if ($ji_shu){
		print "Parameters are not complete.\nWould you like to continue? [Y/N]";
		my $a=<>;
		exit if $a !~ /^YES|Y$/i;
	}
    if ($is_detail){
        if (not $outdir){
            print "OUTDIR has to be given if --detail is specified. Exiting ......" . "\n";
            exit;
        }
        else{
            if (-d $outdir){
                if ($is_force){
                    `rm -rf $outdir`;
                }
            }
            else{
                `mkdir -p $outdir`;
            }

            if ($is_auto_outfile){
                my $suffix = "";
                if ($ka_swi){
                    $suffix = ".ka";
                }
                elsif ($ks_swi){
                    $suffix = ".ks";
                }
                elsif ($kaks_swi){
                    $suffix = ".kaks";
                }
                $outfile1 = join("/", $outdir, basename($seq_list1).$suffix);
                $outfile2 = join("/", $outdir, basename($seq_list2).$suffix);
            }
        }
    }
}


