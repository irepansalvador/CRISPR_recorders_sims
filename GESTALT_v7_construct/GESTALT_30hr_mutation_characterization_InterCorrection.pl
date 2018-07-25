#!/usr/bin/perl -w
use strict;

##################################################################
# This will read a stats file of the GESTALT project and will characterize
# the mutations in the different targets. 
#
# Each stats file contains the data for a biological replicate.
# I analysed the six biological replicates for the "v7 construct 30hr"
# The output files from this script are lists of mutations and their
# frequencies for each target, which then are analysed with the R 
# script "GESTALT_30hr_mutation_characterization_InterCorrection.R"
# provided in the same folder to get the mean saturation and the mutational 
# outcome of each target
# 
# For the mutations it considers separately point target mutations and 
# inter-target mutations. 
# The ultimate objective of this analysis is to obtain experimental
# parameters to perform simulations and quantify the accuracy of the GESTALT
# project.
# (DATA from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81713)
#
# example of input file ($ARGV[0]): "30hr_1_1x.stats"
#
################ Variables             ###########################

my $i;
my $p;
my $s;

my $last;
my %INTER_end;
my %INTER_COUNT;
my %INTER_size;
my %INTER_start;
my %INTER_type;

################ IO ##############################################
my $input = $ARGV[0];
my $output1 = $input; for (1..4) {chop $output1};
$output1 = $output1."_INTER_corrected_results.txt";

my $output2 = $input; for (1..4) {chop $output2};
$output2 = $output2."_INTRA_corrected_results.txt";

my $output3 = $input; for (1..4) {chop $output3};
$output3 = $output3."_OVERALL_corrected_results.txt";

############### Get INTER-TARGET DELETIONS ########################
$i = 0; my $ii = 0;
my $l; my $ll; my $lll;
$last= "NAN";
my $N_inter = 0;

open (FILE1, "$input");
while (<FILE1>)
  {
    chomp;
    my @parts = split "\t", $_;
    #print "$_\n";

    ######### if the line has a valid sequence
  if ($parts[1] =~ m/PASS/)
		{	
			#print "###########\n";
			########## for each target
		for $p (4..13)
			{
			$l = $p - 3;
			#print "$parts[$p]\t";
			my @parts2 = split /&/,$parts[$p] ;
			######### for each mutation (in case there is multiple muts in one target "&")
			for my $p2(0..$#parts2)
				{
				if ($parts2[$p2] =~ m/(\d+)(D|I)/)
					{
					my $n = $1;
					########### if the mutation is a deletion of more than 27pb
					if ($n >= 27)
						{
						if ($N_inter == 0) {$ll = $l; $last = $parts2[$p2]; $ii = $i};
												
						#print "\n read $i: $n FOUND $parts2[$p2] in $parts[$p] in target $l: $N_inter\n";
						
						
						if ($last ne $parts2[$p2] || $ii ne $i)
							{
							if ($N_inter<=1)
								{
								#print "Inter-targets = $N_inter not really! \n";
								### reINTER_start variables
								$ll = $l;
								#$N_inter= 0;
								$last=$parts2[$p2];
								$ii = $i;
								next;
								}
							
							$lll = $ll + $N_inter -1;
							$s = $lll - $ll + 1;
							#print "\n\t found something! from $ll to $lll \n\n";
								### assign the intertarget to the array
							$INTER_COUNT{$last}++;
							$INTER_size{$last}  = "$s";
							$INTER_start{$last} = "$ll";
							$INTER_end{$last} = "$lll";
							$INTER_type{$last} = "$ll-$lll";
							
							#print "INTER of $last : $INTER{$last}\n";
							### reINTER_start variables
							$ll = $l;
							$N_inter= 0;
							$last=$parts2[$p2];
							$ii = $i;
							}
						#print "last : $last N_inter: $N_inter l: $l ll: $ll i : $i  ii : $ii \n";		  
						$N_inter++;
						}
					###########
					}
				}  
			########
			}
		#print "\n";
		$i++;
		#if ($i >= 100) {last}
		}
  }
close FILE1;

print "I have found $i PASS seqs\n";
##########################################################################################
open (OUT1, ">$output1");

print OUT1 "MUTATION\tCount\tStart\tEnd\tSize\tType\n";

foreach my $int (reverse sort { $INTER_COUNT{$a} <=> $INTER_COUNT{$b} } keys %INTER_COUNT)
  {
  print OUT1 "$int\t$INTER_COUNT{$int}\t$INTER_start{$int}\t$INTER_end{$int}\t$INTER_size{$int}\t$INTER_type{$int}\n";
  }

close OUT1;


############### Get INTRA-TARGET DELETIONS ########################
$i = 0; $ii = 0;
$l=0; $ll=0; $lll=0;
$last= "NAN";
$N_inter = 0;

my @NONE =0;
my %INTRA_end;
my %INTRA_COUNT;
my %INTRA_size;
my %INTRA_start;
my %INTRA_type;
open (FILE1, "$input");
while (<FILE1>)
  {
    chomp;
    my @parts = split "\t", $_;
    #print "$_\n";

    ######### if the line has a valid sequence
  if ($parts[1] =~ m/PASS/)
	{	
    #print "###########\n";
    ########## for each target
    for $p (4..13)
      {
      $l = $p - 3;
      #print "$parts[$p]\t";
      my @parts2 = split /&/,$parts[$p] ;
      ######### for each mutation (in case there is multiple muts in one target "&")
	  for my $p2(0..$#parts2)
		{
		if ($parts2[$p2] =~ m/(\d+)(D|I)/)
		  {
		  my $n = $1;
		  $last = $parts2[$p2];
		  ########### if the mutation is a deletion of more than 20pb
		  if ($n <= 28)
			{
			if (exists $INTER_COUNT{$parts2[$p2]})
			  {
				#print "$parts2[$p2] exists in INTER hash\n";
				next;
			  }
			#print "\n read $i: $n FOUND $parts2[$p2] in $parts[$p] in target $l: $N_inter\n";

  			  ### assign the intertarget to the array
			  $INTRA_COUNT{$last}++;
			  $INTRA_start{$last} = "$l";
        
			}
		  ###########
		  }
        elsif ($parts2[$p2] =~ m/NONE/) {$NONE[$l]++}
		}
	  ########
	  }
	#print "\n";
	$i++;
	#if ($i >= 100) {last}
	}
  else {next;}
  }
close FILE1;

#

###########################################################################################
open (OUT2, ">$output2");

print OUT2 "MUTATION\tCount\tTarget\n";

foreach my $int (sort { $INTRA_start{$a} <=> $INTRA_start{$b} } keys %INTRA_start)
  {
  print OUT2 "$int\t$INTRA_COUNT{$int}\t$INTRA_start{$int}\n";
  }

close OUT2;



###########################################################################################
open (OUT3, ">$output3");

print OUT3 "TARGET\tMuts\tReads\tProportion\n";
for (1..10)
    {
    my $prop =  1 - $NONE[$_]/$i;
    print OUT3 "$_\t$NONE[$_]\t$i\t$prop\n";    
    }
print "I have found $i PASS seqs and @NONE NONEs\n";

close OUT3;
