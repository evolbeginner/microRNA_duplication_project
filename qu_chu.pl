#! /usr/bin/perl


$file1=$ARGV[0];
$file2=$ARGV[1];


open(IN,'<',$file1);
while(<IN>){
	chomp;
	@a{(split)} = (1,1);
}
close IN;

open(IN2,'<',$file2);
while(<IN2>){
	chomp;
	@b=split;
	next if exists $a{$b[0]};
	print $b[0]."\t".$b[1]."\n";
}



