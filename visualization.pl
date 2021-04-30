#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;

my $opts = parse_params();
preprocess_data($opts);
generate_html();
generate_js();
create_final_files($opts);

exit;

#--------------------------------

sub error
{
    my (@msg) = @_;
    if ( scalar @msg ) { print "@msg\n"; }
    print 
        "About: This is a script that enables a visualization of results obtained from \"bcftools roh\"\n",
        "Usage: rohViz.pl [OPTIONS]\n",
        "Options:\n",
        "   -i, --input-file <file>       Unzipped input file with output from \"bcftools roh\"\n",
        "   -l, --min-length <num>        Mimimum length of ROH [1e6]\n",
        "   -o, --outdir <dir>            Output directory\n",
        "   -h, -?, --help                This help message\n",
        "\n";
    exit -1;
}

sub parse_params
{
    my $opts =
    {
        min_length  => 1e6,
    };
    while (defined(my $arg=shift(@ARGV)))
    {
		if ( $arg eq '-i' || $arg eq '--input-file' ) { $$opts{in_file}=shift(@ARGV); next }
	    if ( $arg eq '-l' || $arg eq '--min-length' ) { $$opts{min_length}=shift(@ARGV); next }
        if ( $arg eq '-o' || $arg eq '--outdir' ) { $$opts{outdir}=shift(@ARGV); next }
        if ( $arg eq '-?' || $arg eq '-h' || $arg eq '--help' ) { error(); }
        error("Unknown parameter \"$arg\". Run -h for help.\n");
    }
    if ( !exists($$opts{outdir}) ) { error("Missing the -o option.\n"); }
	if ( !exists($$opts{in_file}) ) { error("Missing the -i, --input-file option.\n"); }
	if ( exists($$opts{outdir}) && ! -d "$$opts{outdir}" ) { error("The output directory does not exist: $$opts{outdir}\n"); }
	return $opts;
}

sub preprocess_data
{
	my ($opts) = @_;


	my $dir = "temp";
	if ( !mkdir $dir ) { error("Unable to create temporary directory\n"); }

	for ( my $index = 1; $index < 24; $index++ ) {
		my $name = "temp/chr$index.txt";
		open(OUTPUT, ">$name") or error("Unable to create temporary file\n");
		print OUTPUT "midPoint,sample,snpCount,rateAZ\n";
		close(OUTPUT);
	}

	my $filteredFile = `grep '^[RS]' $$opts{in_file} | sort -n -k3,3 -k2,2 -s > temp/viz$$opts{in_file}`;

	open(ROHOUTPUT, ">temp/roh.txt") or error("Unable to create temporary file\n");

	open(DATA, "<temp/viz$$opts{in_file}") or error("Unable to open the input file\n");

	open(SAMPLES, ">temp/sampleCount.txt") or error("Unable to create temporary file\n");
	
	my $upperBound = 100000;
	my $homozygSNP = 0;
	my $countSNP = 0;
	my $curr_chrom = -1;
	my $curr_sample;
	my $rate;
	my $midPoint;
	my $sampleCount = 0;
	my @splitLine;

	while (my $line = <DATA>) 
	{
		if (index($line,'#') == 0) { next; }
	
		#Split the line
		@splitLine = split(/\s+/, $line);
	
		#ROH
		if ( $splitLine[0] eq "RG" ) 
		{
			if ( int($splitLine[4]) - int($splitLine[3]) >= $$opts{min_length}) 
			{ 
				print ROHOUTPUT "{chr: $splitLine[2], start: $splitLine[3], end: $splitLine[4], sample: \"$splitLine[1]\"},\n"; 
			}
			next;
		}

		#Next chromosome
		if	( int($splitLine[2]) != $curr_chrom )
		{
			if ( $curr_chrom == -1 ) 
			{
				$curr_chrom = int($splitLine[2]);
				$sampleCount = 1;
				open(OUTPUT, ">>temp/chr$curr_chrom.txt");
			} 
			else 
			{
				$curr_chrom = int($splitLine[2]);
				$midPoint = $upperBound - 50000;
				if ( $countSNP == 0 ) { $rate = 0;}
				else { $rate = $homozygSNP / $countSNP; }
				print OUTPUT "$midPoint,$curr_sample,$countSNP,$rate\n";

				close(OUTPUT);

				print SAMPLES "$sampleCount\n";

				$upperBound = 100000;
				$homozygSNP = 0;
				$countSNP = 0;
				$sampleCount = 1;
				open(OUTPUT, ">>temp/chr$curr_chrom.txt");
			}
			$curr_sample = $splitLine[1];
		}

		#Next sample
		if	( $splitLine[1] ne $curr_sample )
		{
			$midPoint = $upperBound - 50000;
			if ( $countSNP == 0 ) { $rate = 0;}
			else { $rate = $homozygSNP / $countSNP; }
			print OUTPUT "$midPoint,$curr_sample,$countSNP,$rate\n";

			$upperBound = 100000;
			$homozygSNP = 0;
			$countSNP = 0;
			$sampleCount++;
			
			$curr_sample = $splitLine[1];
		}

		#Upper bound reached
		if ( int($splitLine[3]) > $upperBound )
		{
			$midPoint = $upperBound - 50000;
			if ( $countSNP == 0 ) { $rate = 0;}
			else { $rate = $homozygSNP / $countSNP; }
			print OUTPUT "$midPoint,$curr_sample,$countSNP,$rate\n";

			$upperBound += 100000;
			$homozygSNP = 0;
			$countSNP = 0;
		}

		#SNP to process
		$countSNP++;
		if ( $splitLine[4] eq "1" ) { $homozygSNP++; }
	}

	#printing last record
	if ( $countSNP == 0 ) { $rate = 0;}
		else { $rate = $homozygSNP / $countSNP; }
	print OUTPUT "$midPoint,$curr_sample,$countSNP,$rate\n";
	print SAMPLES "$sampleCount\n";


	close(OUTPUT);
	close(SAMPLES);
	close(DATA);
	close(ROHOUTPUT);
}

sub generate_html
{
	my $command = `cp vizParts/htmlStart.txt temp/index.html`;
	my $line;
	
	open(HTML, ">>temp/index.html");
	open(SAMPLES, "<temp/sampleCount.txt");

	my $firstFile = 1;
	for ( my $index = 1; $index < 24; $index++ ) {
		my $com2 = `wc -l temp/chr$index.txt`;
		my @splitLine = split(/\s+/, $com2);
		if ( $splitLine[0] > 1 ) { 
			chomp( $line = <SAMPLES> );
			if ( $firstFile )
			{
				print HTML "\t<script src=\"chrViz.js\" type=\"text/javascript\" which_element=\"chrViz$index\" chrom=\"$index\" first=\"Y\" sample_count=\"$line\"></script>\n";
				$firstFile = 0;
			} else {
				print HTML "\t<script src=\"chrViz.js\" type=\"text/javascript\" which_element=\"chrViz$index\" chrom=\"$index\" first=\"N\" sample_count=\"$line\"></script>\n";
			}
		}
	}
	
	print HTML "</body>\n</html>";

	close(SAMPLES);
	close(HTML);
}

sub generate_js
{
	my $command = `cat vizParts/jsStart.txt > temp/chrViz.js`;
	
	open(JS, ">>temp/chrViz.js");

	for ( my $index = 1; $index < 24; $index++ ) {
		my $name = "temp/chr$index.txt";
		my $com2 = `wc -l temp/chr$index.txt`;
		my @splitLine = split(/\s+/, $com2);
		if ( $splitLine[0] == 1 ) { next; }
		open(INPUT, "<$name");
		my $arrIndex = $index - 1; 
		print JS "stringDataArray[$arrIndex] = `";
		while ( my $line = <INPUT>)
		{
			print JS "$line";
		}
		close(INPUT);
		print JS "`;\n\n";
	}

	print JS "var rohs = [";
	open(ROHINPUT, "<temp/roh.txt");
	while ( my $line = <ROHINPUT>)
		{
			print JS "$line";
		}
	close(ROHINPUT);
	print JS "];";
	close(JS);

	$command = `cat vizParts/jsEnd.txt >> temp/chrViz.js`;
}

sub create_final_files
{
	my ($opts) = @_;

	my $command = `cp temp/index.html $$opts{outdir}/index.html`;
	$command = `cp temp/chrViz.js $$opts{outdir}/chrViz.js`;
	$command = `rm -R temp/`;
}