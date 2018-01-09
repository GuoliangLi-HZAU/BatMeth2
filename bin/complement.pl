#!/usr/bin/perl 
#sieves tags into three files
# .des for description, .tag for tags, .phred for score...
# tr -d '\r' <rikky.fq.flat to remove lf..
unless($ARGV[0])
{
	print "Complement bases in sequence data...\n";
	print "tagformat file count.. \n";
	exit;
}

$FileName=$ARGV[0];
unless(open(FILE,$FileName))
{
	print "File does not Exist \n";
	exit;
}

open (TAGFILE,'>'.$FileName.".cmp") or die;
$Line_Count = 0;
$count = 0;
while ($Line = <FILE>)
{
	unless(substr($Line,0,1) eq '>') 
	{
		$Line=~ tr/ACGTacgt/TGCAtgca/;
	}
	print TAGFILE $Line;
}
print "Finished conversion... \n";
close FILE;
close TAGFILE;
