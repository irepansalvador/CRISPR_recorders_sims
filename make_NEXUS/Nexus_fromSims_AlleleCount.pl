#!/usr/bin/perl -w
use strict;
##################################################################
# It takes a simulation from MATLAB and:
# 1) gets some parameters from the simulation file (number of targets, number of states, number of cells)

#Makes a sample of 1000 cells and iterates 100 times:
	# 2) Creates a nexus file to be read by PAUP (creating a character transition matrix)
##################################################################
srand(123); # Initialize the random number generator
################ IO ################################################
my $cells = 0;
my $sample;
my $input = $ARGV[0];
my $weight = $ARGV[1]; #  an argument to specify if the weighted option should be done

print "my input is $input\n";
print "my weight is $weight\n";
my @PATH =  split /\//, $input;
my $sim_file = $PATH[2];
my @parts = split /_/, $PATH[2];
# For creaing NEWICK tree file with the same name just changing the extension 
my $backboneFILE = $sim_file; for (1..4) {chop $backboneFILE};
my $nexus;
my $nexus_dir = "../NEXUS/";	## these folders need to be previously created
my $output_dir = "../OUTPUT/";
my $PAUP_tree;
#Open a file where to store the allele number of each sample
my $Alleles_file = $output_dir.$backboneFILE."_ALLELES";
open (ALLELES, ">$Alleles_file");
close (ALLELES);

open (INPUT, "<$input") or die "cant opppen $input\n"; # the file created with the MATLAB simulation
while (<INPUT>) {$cells++;}
close (INPUT);


## the number of simulation (repetition) useful for output
my $rep = $parts[5]; for (1..3) {chop $rep};

############### defining parameters for analysis ##################
print "\n1) getting parameteres from $ARGV[0]\n";

## if the number of cells is > 1000 just sample 1000 random cells
my $div = $parts[0];
print "my div = $div\n";
($div) = ($div =~ m/(\d+)/);
print "my div = $div\n";
my $number_taxa;
if ($div >=10) {$number_taxa = 1000 +1} #Also will add root
else {$number_taxa= (2**$div) +1}

my $n_samples; # determine the number of cells (max 1K, 10 random samples for div >=10)
if ($div >=10) {$n_samples = 10}
else {$n_samples = 1}

#my $number_taxa = $cells +1;   
my $n_characters = $parts[1];   # Number of targets!
my $n_states = $parts[3] + 1;   	# Number of states!
my $i; my $j;
## define root
my $root = '';
for ($i=1; $i<=$n_characters; ++$i)
	{$root = $root.'0'}
#print "my root = $root\n";

##########  Get the Nmers and create weighted matrix (it can be disabled if the matrix is not specified to be used)  ############
my $freqfile = "../simulations/".$backboneFILE."_freqs.txt";
open (INPUT, "<$freqfile") or die "cant open $freqfile\n"; # the file created with the MATLAB simulation
#my $M = $number_taxa -1; # The number of lines to extract is provided as an argument
my @list = ();
my $contents = "";
my @freq; my @Nmers;
$i = 0;

while (<INPUT>)
	{
	chop $_;
	my @p = split /,/,$_;
	#print "\"$p[1]\",";
	$Nmers[$i] = $p[0];
	$freq[$i]  = $p[1];
	$i++;	
	}
close INPUT;

my $Tot_Nmers = $#Nmers +1;
#foreach(@Nmers) {print "$_ \n"};
##############################################################################
############# SAMPLE 10 TREES WITH 1000 LEAVES #############################
##############################################################################

for $sample (1..$n_samples)
	{
#	print "SAMPLE #$sample\n";
	##################################################################
	if ($weight == 1)
	{
	$PAUP_tree = $backboneFILE."_sample_".$sample."_Weighted_PAUP.tre";
	$nexus = $nexus_dir.$backboneFILE."_sample_".$sample."_Weighted.nexus";		
	}
	if ($weight == 0)
	{
	$PAUP_tree = $backboneFILE."_sample_".$sample."_Unweighted_PAUP.tre";
	$nexus = $nexus_dir.$backboneFILE."_sample_".$sample."_Unweighted.nexus";		
	}

	print "\t 2) creating file = $nexus\n";
	open (NEXUS, "> $nexus");
	print NEXUS "#NEXUS
	
begin data;
	dimensions ntax = $number_taxa  nchar=$n_characters;
	format
	RESPECTCASE
	GAP = -
";
print NEXUS "SYMBOLS = \"@Nmers\"";
print NEXUS "datatype=standard;
matrix

";
	print NEXUS "root\t$root\n";
	#####################################################################
	########### select 100 lines at random ######################
	open (INPUT, "<$ARGV[0]") or die "cant open INPUT\n"; # the file created with the MATLAB simulation
	my $M = $number_taxa -1; # The number of lines to extract is provided as an argument
	my @list = ();
	 
	while (<INPUT>)
		{
		push(@list, $_), next if (@list < $M);
		$list[rand(@list) ] = $_ if (rand($./$M) < 1);
		}
	print NEXUS @list;
	

	close (INPUT);
	print NEXUS ";
end;
";

	##################################################################
	#this is where we print out the end of the nexus file incl assumptions 
	#block and paup block to run tree search
#	###################################################################
#	#############  CALCULATE THE WEIGTHED TABLE     ###################
	## The table should be calculated using the real frequencuies of each
	## Nmer. Each pairwise distance should be 1 over the mean frequency of
	## the Nmer pair.
	if ($weight == 1)
		{
		print NEXUS "begin assumptions;
		[ancstates ancestor = 0:1-$n_characters;]
		usertype mymatrix (stepmatrix)=$Tot_Nmers\n";
		#print("@freq\n");
		print NEXUS "\t@Nmers ";
		print NEXUS "\n";
		
		for $i (0..$#Nmers)
			{
			my $dist=0;my $dist2=0;
			print NEXUS "\t";	
			#print "$#Nmers @Nmers\n";
			for $j (0..$#Nmers)
				{

				if ($j==0)  {$dist = sprintf("%.2f", log(1/($freq[$i])))};
				if ($i==0)  {$dist = sprintf("%.2f", log(1/($freq[$j])))}
				
				if ($j > 0 && $i >0) {$dist = sprintf("%.2f", log(1/($freq[$i])) + log(1/($freq[$j])) )};
					#$dist2= sprintf("%.2f", 2 * $dist);
				if ($j == $i) 			{print NEXUS "0 "}
				#elsif ($j == $n_states)	{print NEXUS "i "}
				elsif ($j==0) {print NEXUS "$dist "}
				else {print NEXUS "$dist "}
				}
				print NEXUS "\n";
			}
			#}
	
		print NEXUS ";
	typeset mytypes = mymatrix: 1-$n_characters;
end;
";
	
	print NEXUS "begin paup;
	ctype mymatrix: all;
	nj brlens = yes;
	savetrees format = altnexus file = $PAUP_tree replace = yes;
	quit nowarn;
end;
";

		}
	##---------------------------------------------------------------------------		
	if ($weight == 0)
	{
print NEXUS "begin paup;
	nj brlens = yes;
	savetrees format = altnexus file = $PAUP_tree replace = yes;
	quit nowarn;
end;";		
	}

	##################################################################
	# Count the number of alleles from the recent NEXUS file
	##################################################################

	open (NEXUS, "<$nexus") or die "cant open NEXUS\n";
	my $cont  =0; 	my $line;
	my @unique;
	# open the file with 1000 sampled cells of the GESTALT
	while (<NEXUS>)
    {
    chop $_;    $line=$_;
		# get the lines with the cells info
    if ($line =~ m/c_00/)   
    	{
			my $found=0;
			## get the mutated targets
    	my @parts = split /\t/, $line;
    	#--
			foreach (@unique) {if ($parts[1]=~$_)	{$found=1}}
			## if is not found, add the mutated target to the unique array
			if ($found==0) {push @unique,$parts[1];
											$cont++}
	    }
	  }
	close (NEXUS);
	## print the results
	#foreach (@unique) {print "$_\n"}
	open (ALLELES, ">>$Alleles_file");
	print ALLELES "$cont\n";
	close (ALLELES);
	}	
