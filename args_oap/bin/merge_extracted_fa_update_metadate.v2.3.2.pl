#!/usr/bin/perl -w
use strict;

#Authur: Xiao-Tao Jiang
#Email biofuture.jiang@gmail.com

##Generate input for online part analysis, merge all extracted fasta 
##and update the meta_data file with number of reads for each sample and number of 16S like reads in each sample

die "perl $0 <Indir> <outdir> <meta_data_in> <meta_data_out> <extracted fasta>  <option_c> <coglist> $!\n" unless (@ARGV == 7);

use Getopt::Std;
use File::Basename;
use FindBin qw($Bin);

##Generalize dir for this program
our (@dirset,$ublastxdir);
BEGIN {
    @dirset = split(/\//,$Bin);
    $ublastxdir = join("/", @dirset);
    unshift @INC, "$ublastxdir/bin";
}


#Indir is the directory including your fastq file
#Outdir is the directory for output of ublastx_v1.1
#meta_data_in is the orignial meta_data file
#meta_data_out is the update meta data file for online ARG analysis
#extracted fasta is the merged extracted potential ARGs reads for online analysis

##option_c is S/U
##coglist is  all_KO30_name.list

#copycorrect is the indicator for whehter perform copy number correction


my %sampleid;
my %metainfo;
my %samplerlen;

#my $coglis = "../DB/all_KO30_name.list";
my $coglis = $ARGV[6];

die "$!\n" unless open(I,"$ARGV[2]");
die "$!\n" unless open(T,">$ARGV[3]");
die "$!\n" unless open(TT,">$ARGV[4]");

my $h = <I>; chomp($h);
print T "$h\t#ofReads\t#of16Sreads\tCellNumber\n";

while(<I>){
	chomp;
	my @tm = split("\t", $_);
	$sampleid{$tm[0]} = $tm[1];
        $samplerlen{$tm[0]} = $tm[3];
	$metainfo{$tm[0]} = $_;
}
close I;


my $opt_i = $ARGV[0];
my $opt_o = $ARGV[1];
my %hashreads;
my %hash16s;
my %cellnum;

for my $ids (sort {$a <=> $b}  keys %sampleid){

	##for each sample statistic the number of reads and number of 16S reads
	if(-e "$opt_i/$sampleid{$ids}_1.fa" && -e "$opt_i/$sampleid{$ids}_1.fa"){
			my $num1 = `grep ">" $opt_i/$sampleid{$ids}_1.fa -c`;
			my $num2 = `grep ">" $opt_i/$sampleid{$ids}_2.fa -c`;
			my $ns = $num1 + $num2;
			$hashreads{$sampleid{$ids}} = $ns;
	}else{
		die "Can not find the input files\n";
	}
	
	##for 16s
	my $s16num1s = 0;
	if(-e "$opt_o/$sampleid{$ids}_1.sam"){
               my $s16num1s = `$ublastxdir/samtools view -F 4 $opt_o/$sampleid{$ids}_1.sam|wc -l`;
        }else{
               die "Can not find the output of _1.fa from minimap2\n";
        }

       	if(-e "$opt_o/$sampleid{$ids}_2.sam"){
               my $s16num2s = `$ublastxdir/samtools view -F 4 $opt_o/$sampleid{$ids}_2.sam|wc -l`;
	       my $s16nums = $s16num1s + $s16num2s;
               $hash16s{$sampleid{$ids}} = $s16nums * $samplerlen{$ids} / 1432;
        }else{
	       $hash16s{$sampleid{$ids}} = $s16num1s * $samplerlen{$ids} / 1432;
               die "Can not find the output of _2.fa from minimap2\n";
        }


        ##for cell number counting from the USCMGs 
	if(-e "$opt_o/$sampleid{$ids}.uscmg.blastx.txt"){
        my $avgKOcopy = 0;        
        my %seq2OGs;
        my %seqlen;
        die "$!\n" unless open(OGMAP, "$coglis");
        while(<OGMAP>){
            chomp;
            my @tem = split(/\t/, $_);
            $seq2OGs{$tem[0]} = $tem[1];
            $seqlen{$tem[0]} = $tem[2];
        }
        close OGMAP;
        die "$!\n" unless open(TAB, "$opt_o/$sampleid{$ids}.uscmg.blastx.txt");
        die "$!\n" unless open(KOAVERAGE, ">$opt_o/$sampleid{$ids}.uscmg.ko_averagecov.txt");
        my %seqcov;
        while(<TAB>){
            chomp;
            my @tem = split(/\t/, $_);
            if(exists $seq2OGs{$tem[1]}){
                $seqcov{$tem[1]} += $tem[3];
            }else{
                $seqcov{$tem[1]} = $tem[3];
            }
			}
        close TAB;            
        my %kocov;
        for my $sid (keys %seqcov){
            if(exists $seq2OGs{$sid}){
                if(exists $kocov{$seq2OGs{$sid}}{$sid}){
                    $kocov{$seq2OGs{$sid}}{$sid} += $seqcov{$sid};
                }else{
                    $kocov{$seq2OGs{$sid}}{$sid} = $seqcov{$sid};
                }
            }else{
                die "Wrong ID $sid\n";
            }
        }
        my $KOcount = 0;
        for my $koid (keys %kocov){
            my $ave = 0; 
            my $seqnum =0;
            for my $seqid (keys %{$kocov{$koid}}){
                my $seqaverage =  $kocov{$koid}{$seqid} / $seqlen{$seqid};
                $ave += $seqaverage;
                $seqnum ++;
            }
            #$ave = $ave/ $seqnum;
            print KOAVERAGE "$koid\t$ave\t$seqnum\n";
            $avgKOcopy += $ave;
            $KOcount ++;
        }
        die "Line 154 $0 Not KO Maping, wrong datasets\n" if($KOcount == 0);				
        $cellnum{$sampleid{$ids}} = $avgKOcopy / $KOcount; ##The aveage coveage of each universal single copy marker genes KO is the estimated cell number 
    }else{
    die "Wrong run, without the merge blastx result\n" 
    }
	##
    	##merge pair_end file reads into one with the sample IDs
    	if(-e "$opt_o/$sampleid{$ids}.extract_1.fa" && -e "$opt_o/$sampleid{$ids}.extract_2.fa"){
    		die "$!" unless open(TEM1, "$opt_o/$sampleid{$ids}.extract_1.fa");
    		die "$!" unless open(TEM2, "$opt_o/$sampleid{$ids}.extract_2.fa");
    
    		my $count = 1;
    		while(my $ntem = <TEM1>){
    			my $stem = <TEM1>;
    			my $changename = ">$sampleid{$ids}_$count";
    			print TT "$changename\n$stem";
    			$count ++;
    		}
    		close TEM1;
    	
    		while(my $ntem = <TEM2>){
    			my $stem = <TEM2>;
    			my $changename = ">$sampleid{$ids}_$count";
    			print TT "$changename\n$stem";
    			$count ++;
    		}
    		close TEM2;
    		$count = 0;
    
    	}else{
    		die "Can not find the extracted.fa\n";
    	}		

}##process each pair of fastq 

for my $id (sort {$a <=> $b} keys %metainfo){
	##update meta data files 
	print T "$metainfo{$id}\t$hashreads{$sampleid{$id}}\t$hash16s{$sampleid{$id}}\t$cellnum{$sampleid{$id}}\n";
}

__END__
1;
