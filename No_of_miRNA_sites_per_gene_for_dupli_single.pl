#! /bin/env perl

use strict;
use Getopt::Long;

#####################################################################
sub get_genes_from_list{
    my %genes;
    my ($gene_list_file) = @_;
    open (my $IN, '<', $gene_list_file) || die "gene_list_file $gene_list_file cannot be opened!";
    while(my $line = <$IN>){
        chomp($line);
        $genes{$line}="";
    }
    return(\%genes)
}

sub get_num_of_target_sites{
    my %num_of_target_sites;
    my ($targets_file) = @_;
    open(my $IN, '<', $targets_file) || die "targets file cannot be opened!";
    while(my $line = <$IN>){
        chomp($line);
        my @line_array = split(/\t/, $line);
        @line_array = grep {$_ ne ""} @line_array;
        $num_of_target_sites{$line_array[0]} = $#line_array;
    }
    return(\%num_of_target_sites);
}

#####################################################################
my ($duplicates_file, $singletons_file, @targets_files);
my (%genes);

GetOptions(
    'd|dupli|duplicates=s'    =>  \$duplicates_file,
    's|single|singletons=s'   =>  \$singletons_file,
    't|target_file=s'         =>  \@targets_files,
) || die "Params error!";

#####################################################################
$genes{'duplicates'} = &get_genes_from_list($duplicates_file);
$genes{'singletons'} = &get_genes_from_list($singletons_file);

for my $targets_file (@targets_files){
    print $targets_file."\n";
    my (%final_results, %counter);
    my $num_of_target_sites_for_each_gene_href = &get_num_of_target_sites($targets_file);
    while (my ($type, $genes_href) = each %genes){
        foreach my $gene (keys %$genes_href){
            if (exists $num_of_target_sites_for_each_gene_href->{$gene}){
                $final_results{$type} += $num_of_target_sites_for_each_gene_href->{$gene};
                $counter{$type} ++;
            }
        }
    }
    while(my ($type, $final_number) = each %final_results){
        print $type."\t".$final_number/$counter{$type}."\n";
    }
    print "\n";
}


