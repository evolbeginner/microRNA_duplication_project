use Statistics::Basic qw(:all);

open(IN,"./宽松尺度的paralog.txt");
while(<IN>)
{
       chomp;
       my @a=split;
       $target{$a[0]}=$a[1];
       $target{$a[1]}=$a[0];
}
close IN;


open(IN,"../../haha2");
#mayer TD 原始文件.txt
while(<IN>)
{
       chomp;
       $hash{$_}=1;
       my @a=split;
       do {$_=uc($_); $element{$_}=1} foreach (@a);
       next if not exists 'target' -> {$a[0]};
       $para{$a[0]}=$a[1];
       $para{$a[1]}=$a[0];
       my $pair=join "-", sort ($a[0], $a[1]);
       $pair{$pair}=1;

}
close IN;


open (IN,"./StrippedMicroData.txt");
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
        print $symbol."\t".$para."\t".correlation(\@x,\@y)."\n";
}









