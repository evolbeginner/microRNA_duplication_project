#! /usr/bin/perl
use 5.010;
use Bio::SeqIO;

my ($filtered_seq,$DNA_file,$ref_seq,$ref_filtered_seq,$pep_file,$ref_longest);
my (%filtered_seq,%ref_filtered_seq);

$DNA_file = "../sequence/TAIR9_cds.fasta";
$pep_file = "../sequence/TAIR9_pep.fasta";


##################################################################
$ref_seq = &read_fasta($DNA_file);

$ref_filtered_seq = &filtered_seq($ref_seq);

{
	$ref_seq1 = &extract_seq($DNA_file,$ref_filtered_seq);
	$ref_seq2 = &extract_seq($pep_file,$ref_filtered_seq);	
}

$ref_longest = &generate_longest_seq_name($ref_filtered_seq);

{
	&output($ref_longest, $ref_seq1, 'TAIR9_DNA_longest.fasta');
	&output($ref_longest, $ref_seq2, 'TAIR9_pep_longest.fasta');
}

##########################################################################
##########################################################################
sub read_fasta{
	my ($input_seq) = @_;
	my ($seq,$seq_obj,$catchseq_seqio_obj);
	my (%seq,$length);
	
	$catchseq_seqio_obj = Bio::SeqIO->new(-file=>$input_seq, -format=>'fasta');
	while($seq_obj=$catchseq_seqio_obj->next_seq()){
		my ($name,$seq);
		$name = $seq_obj->display_name;
		my ($new_name) = $name =~ /\>([^ ]+)/;
		$seq = $seq_obj->seq;
		$length = $seq_obj->length;
		$seq{$name}{'seq'}=	$seq;
		$seq{$name}{'length'}=	$length;
	}
	return (\%seq);
}

sub filtered_seq
{
my ($ref_seq) = @_;
foreach my $seq_name(keys %$ref_seq){
	my ($core_seq_name);
	($core_seq_name) = $seq_name =~ /([^\.]+)\.\d+/; # no \.
	if (not exists $filtered_seq{$core_seq_name}){
		&generate_filtered_seq($core_seq_name,$seq_name);
	}
	else{
		$ref_filtered_seq->{$core_seq_name}{'length'} < $ref_seq->{$seq_name}{'length'} ? &generate_filtered_seq($core_seq_name,$seq_name) : 1 ;
	}
}
return($ref_filtered_seq);
}


sub generate_filtered_seq{
	my ($core_seq_name,$seq_name) = @_;
	$ref_filtered_seq->{$core_seq_name}=$ref_seq->{$seq_name};
	$ref_filtered_seq->{$core_seq_name}{'name'}=$seq_name;
}


sub extract_seq{
	my ($pep_file,$ref_filtered_seq) = @_;
	my (%ref_seq2);
	$ref_seq2=&read_fasta($pep_file);
	return($ref_seq2);
}


sub output{
	my ($ref_longest, $ref, $outfile_name) = @_;
	open(my $out,'>',$outfile_name);
	select $out;
	foreach (sort keys %$ref){
		my $name=$ref->{$_}{'name'};
		next if not exists $ref_longest->{$_};
		print '>'.$_."\n";
		print $ref->{$_}{'seq'}."\n";
	}
	close $out;
}


sub generate_longest_seq_name{
	my ($ref) = @_;
	my (%ref_longest,$core_seq_name);
	foreach my $core_seq_name(keys %$ref){
		$ref_longest -> {$ref->{$core_seq_name}{'name'}} = 1;
	}
	return($ref_longest);
}


