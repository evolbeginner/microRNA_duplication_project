#! /usr/bin/perl

$target = "target_2-100.list";

$dup = $ARGV[0];;
$single = $ARGV[1];

open(IN,'<',$target);
%target = map {chomp;$_ => 1} <IN>;

#foreach(keys %target){print $_."\n";}


foreach($dup,$single)
{
	my $count;
	open(IN,$_);
	while(<IN>){
		chomp;
		$count++ if exists $target{$_};
	}
	print $count."\n";
}



