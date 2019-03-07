#!/usr/bin/perl 
#strip will make a flat file out of a fasta file....

unless($ARGV[0] and $ARGV[1])
{
	print "Strip will create a flat file out of filename 1 and save it as filename 2 \n";
	exit;
}
$FileName=$ARGV[0];
$OutFile=$ARGV[1];


unless(open(FILE,$FileName))
{
	print "File does not Exist \n";
	exit;
}
open (FLATFILE,'>'.$OutFile) or die;
while($Line=<FILE>)
{
	unless(substr($Line,0,1) eq '>') 
	{
		 chomp($Line);
		 print FLATFILE lc($Line);
		 #print substr($Reverse,1);
	}
}
close FILE;
