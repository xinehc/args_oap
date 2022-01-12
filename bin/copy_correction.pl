#!/usr/bin/perl -w
use strict;

die "perl $0 <community> <copy_writer_db> <average_copynum>\n" unless (@ARGV == 3);

die "$!\n" unless open(I, "$ARGV[0]");
die "$!\n" unless open(II, "$ARGV[1]");
die "$!\n" unless open(T, ">$ARGV[2]");


my %copyn;
<II>;
while(<II>){
	chomp;
	my @tem = split(/\t/, $_);
	$copyn{$tem[0]} = $tem[1];
}
close II;

my $total = 0;
my %weight;

while(<I>){
	chomp;
	my @ts = split(/\t/, $_);
	$total += $ts[1];
	if(exists $copyn{$ts[0]}){
	
		my $acnum = $ts[1] * $copyn{$ts[0]};
		#print T "$ts[0]\t$acnum\n"; 
		$weight{$ts[0]} = $acnum;
	
	}else{
		my @tttmp = split("; ", $ts[0]);
		my $val = copyupper(@tttmp);
		#my $acnum = $ts[1] * 2.45;
		my $acnum = $ts[1] * $val;
		#print "$ts[0]\n";
		#print T "$ts[0]\t$acnum\n";
		next if($val == 0);
		$weight{$ts[0]} = $acnum;
	}
}
close I;

my $normacopy =0;
for my $k (keys %weight){

	$normacopy += ($weight{$k} / $total);
}

print T "Normalizedcopynumber\t$normacopy\n";


sub copyupper {
	my @all = @_;
	my @tem = ("k__", "p__", "c__", "o__", "f__", "g__", "s__");
	my $cop = 1;

	for(my $i = 0; $i <= 6; $i++){
		$tem[$i] = $all[$i];
		my $test = join("; ", @tem);
		if(exists $copyn{$test}){
			$cop = $copyn{$test};
		}else{
			return 0;
		}
	}	
	#print "$cop\n";
	return $cop;
}
__END__
1;
