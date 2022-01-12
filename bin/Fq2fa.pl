#!/usr/bin/perl -w
use strict; 
die "perl $0 <Ifq> <Ofa>\n" unless (@ARGV ==2);
die "$!\n" unless open(I, "$ARGV[0]");
die "$!\n" unless open(T, ">$ARGV[1]");

while(my $n=<I>){
	my $se = <I>;
	<I>;<I>;
	$n =~ s/^@/>/;
	print T "$n$se";
}
close I;
close T;
