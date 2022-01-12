#!/usr/bin/perl -w 
use strict;

#Author Xiao-Tao JIANG
#Date 2018/1/6
##------------------------------------
#This script is written to split too
#big extracted.fa into smaller ones
#in order to run ARGs-OAP stage two paralelly
##-----------------------------------

die "perl $0 <extracted.fa> <meta_data_online.txt> <samples_per_part>\n\tThe third parameter is the lines per divided part included, \n\twe suggested for one part the maxmum size is upto 500 Mbs, so if the size of your single extracted.fa is 50 M, about ten samples each part is recommended\n\trun the program like this\n\tExample: perl $0 extracted.fa meta_data_online.txt 10\n" unless (@ARGV ==3);

##Parse the meta_data_online.txt to decide how many parts should be generated. 
my %mot; 
die "$!\n" unless open(MO, "$ARGV[1]");
my $lines = `wc -l $ARGV[1]`;
$lines = (split(/\s+/, $lines))[0] - 1;
my $head = <MO>;
die "Incorrect number of lines per part: $ARGV[2]\n" if($ARGV[2] <= 0);
my $nparts = int($lines / $ARGV[2]);
my $mod = $lines % $ARGV[2];
if($mod == 0){
}else{
	$nparts ++;
}
#print "$nparts\n";

##split the meta_data_online.txt into #nparts
my @samples = ();
while(<MO>){
	chomp;
	push @samples, $_;
}
close MO;

##generate the divided meta data online files 
if($mod == 0){ ##JUST divided OK
	for(my $i = 1; $i <= $nparts; $i++){
		my $ots = "meta_data_online_part$i.txt";
		die "$!\n" unless open(OTN, ">$ots");
		print OTN "$head";
		for(my $j = ($i-1) * $ARGV[2]; $j < $i * $ARGV[2]; $j++){
			print OTN "$samples[$j]\n"; 
		}	
		close OTN;
	}
}else{
	##The last one does not have all the number of samples 
	if($nparts == 1){
			my $ots = "meta_data_online_part1.txt";
			die "$!\n" unless open(OTN, ">$ots");
			print OTN "$head";
			for(my $j = 0;  $j <= $#samples; $j++){
				print OTN "$samples[$j]\n"; 
			}	
			close OTN;

	}elsif($nparts > 1){
		#first $nparts - 1 parts 
		for(my $i = 1; $i < $nparts; $i++){
			my $ots = "meta_data_online_part$i.txt";
			die "$!\n" unless open(OTN, ">$ots");
			print OTN "$head";
			for(my $j = ($i-1) * $ARGV[2]; $j < $i * $ARGV[2]; $j++){
				print OTN "$samples[$j]\n"; 
			}	
			close OTN;
		}
		##Last part output 
		my $ots = "meta_data_online_part$nparts.txt";
                die "$!\n" unless open(OTN, ">$ots");
		print OTN "$head";
		for(my $j = ($nparts -1) * $ARGV[2]; $j <= $#samples; $j++){
				print OTN "$samples[$j]\n";
		}
		close OTN;
	}else{

		die "Wrong Number $nparts\n";
	}
}#OUTPUT THE DIVIDED META_DATA_ONLINE.TXT
#the last one

##Check whether the dividing process is correct or not
for(my $i=1; $i<=$nparts; $i++){
	die "Without Parts $i\n" unless(-e "meta_data_online_part$i.txt");
}

##Generate the divided extracted.fa for each part

for(my $i =1; $i<=$nparts; $i++){
	
	my $otfa = "extracted_fa_part$i.fa";
        die "$!\n" unless open(PFA, ">$otfa");
	
	my %namehash;
	die "$!\n" unless open(MA, "meta_data_online_part$i.txt");
	<MA>;
	while(<MA>){
		chomp;
		my @m = split /\t/;
		$namehash{$m[1]} = 1;
	}
	close MA;

	die "$!\n" unless open(FA, "$ARGV[0]");
	while(my $n = <FA>){
		my $s = <FA>;
		my $bn = $n;
		chomp($n); $n =~ s/^>//;
		my @m = split(/\s+/, $n);
		
		$m[0] =~ s/\_(\d+)$//; ##fecth the sample name 

		if(exists $namehash{$m[0]}){
			print PFA "$bn$s";
		}
	}
	close FA;
	close PFA;
	##delete %namehash
	%namehash = ();
}

##FINISH DIVIDING
