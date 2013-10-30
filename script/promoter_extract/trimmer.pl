#!/usr/local/bin/perl

#open SORTEDINPUT, "sorted.bed" or die "Cannot open sorted file.\n";

@lines = <STDIN>;

chomp(@lines);

$filesize = @lines;

#close SORTEDINPUT;

#open OUTPUTFILE, ">query_results.bed" or die "Could not create output file.\n";

foreach $line (@lines)
{
	@temp = split('\t', $line);

	$size = @temp;

	for($i = 0; $i < $size; $i++)
	{
		if($i != 1)
		{
			$tempprint = $temp[$i];
			print $tempprint;
			if(($i+1) < $size)
			{
				print "\t";
			}
		}
	}
	print "\n";	
}

#close OUTPUTFILE;
