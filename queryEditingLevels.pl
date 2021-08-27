#!/bin/perl -w
############################################################
#Original Author: Gokul 
#perl script that queries editing level of known sites in a BAM list file

use strict;

my ($inputfile, $outputfile) = ($ARGV[0], $ARGV[1]);

#GLOBAL VARIABLES - PLEASE MODIFY THESE ACCORDINGLY

my $minbasequal = 30; # MINIMUM BASE QUALITY SCORE
my $minmapqual = 0; # MINIMUM READ MAPPING QUALITY SCORE
my $sampath = "samtools"; #PATH TO THE SAMTOOLS EXECUTABLE
my $offset = 32; #BASE QUALITY SCORE OFFSET - 32 FOR SANGER SCALE, 64 FOR NEW ILLUMINA SCALE
my $reference = "/ref/autosomes.hg38.fa" # REFRENCE FATSTA FILE LOCATION

##END GLOBAL VARIABLES

open (my $INPUT , "<", $inputfile) or die "error opening inputfile: $!\n"; #INPUT AN RNA EDITING REFERNCE PANEL 
open (my $OUTPUT1, ">", "$outputfile\.cov"); # COVERAGE READS NUMBER
open (my $OUTPUT2, ">", "$outputfile\.edi"); # EDITING READS NUMBER
#open (my $OUTPUT3, ">", "$outputfile\.edfreq"); 
#print $OUTPUT "#chrom\tposition\tinfo\tstrand\tcoverage\teditedreads\teditlevel\n";

my %pos2info=();
my %pos2strand=();
while (<$INPUT>) { #READ IN LIST OF KNOWN EDITED SITES AND QUERY EDITING STATUS
$_=~s/\s+$//;
my @fields = split(/\t/,$_);
my $pos="$fields[0]\t$fields[2]";
my @siteinfo = split(/\,/, $fields[3]);
my $strand = $siteinfo[1];
$pos2info{$pos}="$pos\t$fields[3]";
$pos2strand{$pos}=$strand;
}
close $INPUT;

#my $TEMPNAME = join '', $outputfile,'_mpileup';

# OPEN A PIPE FORM SAMTOOLS MPILEUP
open(MPILEUP, "$sampath mpileup -x -A -B -C 0 -l $inputfile -d 100000000 -q $minmapqual -Q 0 -R --ff 4096 -f $reference -b $outputfile |");

#open(TMPFILE, $ARGV[1])or die;
while (<MPILEUP>) { #READ MPILEUP FROM STDIN
$_=~s/\s+$//;

my @line=split(/\t/,$_);
my ($chr, $pos, $ref)=@line[0..2];

my $j=3;
my $temp_cov="";
my $temp_mismatch="";
my $temp_varfreq="";
my $strand=$pos2strand{"$chr\t$pos"};
my $info = $pos2info{"$chr\t$pos"};
	if($strand ne ""){
	my $editnuc = 'G';
		if ($strand eq '-'){$editnuc = 'C';}
		
		while(($j+2)<@line){
		my ($counts, $alignment, $qualities) = ($line[$j],$line[$j+1],$line[$j+2]);
		$alignment=~s/\$|\^.//g; #substituted head and tail marker
			while($alignment=~/[+-]([0-9]+)[ACGTNacgtn]+/){
			my $number=$1;
			my $front=$`; #before matched part
			my $match=$&; #matched part
			my $behind=$'; #behind matched part
			$match=~s/[+-]$number[ACGTNacgtn]{$number}//;
			$alignment=$front.$match.$behind;
			}
		my @qualscores = split(//,$qualities);
		my @sequencebases = split(//,$alignment);
		my $i = 0;
		my ($newcov, $newmismatch, $varfreq) = (0,0,"NA");
			while ($i < @qualscores) {
			my $base=$sequencebases[$i];
			my $baseQ=$qualscores[$i];
				if(ord($baseQ)>=($minbasequal+$offset)){ # baseQ >= cutoff
					if($base=~/[\,\.ATCGatcg]/){ # base != N
					$newcov++;
						if ($base=~/$editnuc/i){$newmismatch++;}
					}
				}
			$i++;
			}
			
			if ($newcov!=0) {
			$varfreq = ($newmismatch/$newcov);
			}
			
			if($temp_cov eq ""){
			($temp_cov,$temp_mismatch,$temp_varfreq)=($newcov,$newmismatch,$varfreq);
			}
			else{
			$temp_cov="$temp_cov\t$newcov";
			$temp_mismatch="$temp_mismatch\t$newmismatch";
			$temp_varfreq="$temp_varfreq\t$varfreq";
			}
		$j=$j+3;
		}
	print $OUTPUT1 "$info\t$strand\t$temp_cov\n";
	print $OUTPUT2 "$info\t$strand\t$temp_mismatch\n";
	#print $OUTPUT3 "$info\t$strand\t$temp_varfreq\n";
	}
}
close MPILEUP;
close $OUTPUT1;
close $OUTPUT2;
#close $OUTPUT3; # will calulate editing frequencies later
