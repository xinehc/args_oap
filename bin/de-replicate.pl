#!/usr/bin/perl -w
use strict;

die "perl $0 <Fasta> <No-replicate-annotation> <With-replicate-annotation>\n" unless (@ARGV == 3);
##This script is used to dereplicate ARG sequences with 100% identity (exactly the same sequences)
##and output all the annotation information for that sequence
##This is for the convenient of selection of the most explainable annotation for that ARG sequence

die "$!\n" unless open(I,"$ARGV[0]");
die "$!\n" unless open(T, ">$ARGV[1]");
die "$!\n" unless open(TT, ">$ARGV[2]");

my %unique;
my %rep;
while(my $name = <I>){
        my $seq = <I>;
        if(exists $rep{$seq}{$name}){
                die "Totally replicate $name\n";
        }else{
                $rep{$seq}{$name} = 1;
        }
}
close I;


for my $s (keys %rep){
        if(scalar keys %{$rep{$s}}  == 1 ){
                for my $n (keys %{$rep{$s}}){
                        print T "$n";
                }
                print T "$s";

        }else{

                for my $n (keys %{$rep{$s}}){
                        print TT "$n";
                }
                print TT "$s";

        }
}
