#!/usr/bin/perl -w
use strict;


##This script is to obtain the microbial community from the searched hyper variable sequences


die "perl $0 <hyper.out.fa> <Taxonomy.txt> <community.txt>\n" unless(@ARGV == 3);
die "$!" unless open(I, "$ARGV[0]");
die "$!" unless open(II, "$ARGV[1]");
die "$!" unless open(T, ">$ARGV[2]");

my %id2tax;
while( <II>){
	chomp; 	
	my @tem = split("\t", $_);
	$id2tax{$tem[0]} = $tem[1];
}
close II;

my @blank = ("k__", "p__", "c__", "o__", "f__", "g__", "s__");
sub lca{
        my @all = @_;
        my @total = ();
        my @common = ();

        ##split each taxonomy by "; " 


        my $i = 0;
        my $j = 0;
        ##store in matrix and transform column to row
        for($i = 0; $i <= $#all; $i++){
                my @tem = split("; ", $all[$i]);
                die "$all[$i]" if($#tem != 6);
                for($j = 0; $j <= $#tem; $j++){
                        $total[$j][$i] = $tem[$j];
                }
        }

        ##compare and generate consistent taxonomy assignment by over 2/3 of all the matched
        for(my $k = 0; $k < $j; $k++){
                my $flag = 0;

		my %temp = ();
		
                for(my $l = 0; $l < $i; $l++){
			if(exists $temp{$total[$k][$l]}){
				$temp{$total[$k][$l]} ++;
			}else{
				$temp{$total[$k][$l]} = 1;
			}	
                }

		my @temn = sort {$temp{$b} <=> $temp{$a}} keys %temp;
		if($temp{$temn[0]} >= 2*($k+1)/3){
			$common[$k] = $temn[0];

		}else{
			$common[$k] = $blank[$k];
		}
        }

	return join("; ", @common);
}


my %community;
my $totalv6 = 0;
while(my $n = <I>){
	my $s = <I>;
	chomp($n);
	my @tem = split(/\_/, $n);
	#print "$tem[1]\n";

	my @alls = split("-", $tem[-1]); ##split all the hits to one reads
	my @seqtax;
	my $averlen = 0;
	#my $maxle = 0;
	for(my $i =0; $i <= $#alls ; $i++){
		my @sp = split(":", $alls[$i]);
		$averlen += $sp[1];
		#if($sp[1] >= $maxle){
		#	$maxle = $sp[1];
		#}
		#print "$sp[0]\t$id2tax{$sp[0]}\n";
		push @seqtax, $id2tax{$sp[0]};
	}
	##For abundance
	$averlen = $averlen / ($#alls + 1);
	
	##for taxonomy LCA
	my $taxs = lca(@seqtax);
	#print "$taxs\n";

	#if(exists $id2tax{$tem[-1]})
	if(exists $community{$taxs}){
		$community{$taxs} += $averlen/66.2;
		#$community{$taxs} += $maxle/66.2;
	}else{
		$community{$taxs} = $averlen/66.2;
	}

	#$totalv6 += $maxle; 
}
close I;

#print "$totalv6\n";
for my $key (sort keys %community){
	print T "$key\t$community{$key}\n";
}



