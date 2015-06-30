#! /usr/bin/perl

use Statistics::Basic qw(:all);
use List::Util qw(sum first);
use Statistics::TTest;
use Getopt::Long;

##############################

GetOptions(
	'pair_list_file=s'	=>	\$pair_list_file,
	'pair_seperator=s'	=>	\$pair_seperator,
	'microarray_data=s'	=>	\$microarray_data,
) || die "illegal param";
print $pair_list_file."\n";

die "pair_list_file has to be specified by --pair_list_file" if not $pair_list_file;
$pair_seperator="\t" if not $pair_seperator;
$microarray_data="/home/sswang/project/microRNA/official/data/StrippedMicroData.txt" if not $microarray_data;
#####################################################################################

open(IN,"$pair_list_file") or die "input file cannot be opened";
while(<IN>)
{
	chomp;
	$hash{$_}=1;
	my @line = split/\t/;
	my @a=split /$pair_seperator/, $line[0];
	do {$_=uc($_); $element{$_}=1} foreach (@a);
	$para{$a[0]}=$a[1];
	$para{$a[1]}=$a[0];
	my $pair=join "-", sort ($a[0], $a[1]);
	$pair{$pair}=1;
}
close IN;

open (IN,"$microarray_data") or die "data file cannot be opened";
while(<IN>)
{
       chomp;
       my @a=split;
       $a[0]=uc($a[0]);
       next if not exists $element{$a[0]};
       @{$a[0]}=@a;
       shift @{$a[0]};
}


foreach my $pair(keys %pair)
{
        #next if exists $Cun_Zai{$symbol};
        my ($seq1,$seq2) = split /\-/,$pair;
        @x=@{$seq1};
        @y=@{$seq2};
        next if $#x != $#y;
	next if scalar @x == 0;
	next if (not $seq1) or (not $seq2);
	compare_paralog_pair(\@x,\@y);
        print $seq1."\t".$seq2."\t".correlation(\@x,\@y)."\n";
        #print correlation(\@x,\@y)."\n";
}


sub compare_paralog_pair
{
my ($array_ref1,$array_ref2) = @_;
foreach(@_){
        $mean = &mean(@$_);
        #print $mean."\t";
}
$ttest = new Statistics::TTest;
$ttest->load_data(@_);
my $prob = $ttest->{t_prob};
#print $prob."\n";
}


