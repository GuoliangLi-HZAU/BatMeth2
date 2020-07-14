#!/usr/bin/perl 
print "Converting bad bases to ACGT...\n";

$FileName=$ARGV[0];
unless(open(FILE,$FileName))
{
        print "File does not Exist \n";
        exit;
}

if($ARGV[1]) {open (FILTER,'>'.$ARGV[1]) or die;}
else {open (FILTER,'>'.$FileName.".filter") or die;}
$Line_Count = 0;
$count = 0;
while ($TagLine = <FILE> )
{
        if (substr($TagLine,0,1) eq '>')
        {
                print FILTER $TagLine;
        }
        else
        {
                chomp $TagLine;
                $TagLine=~ tr/acgt/ACGT/;
                $TagLine=~ tr/ACGT/C/c;
                print FILTER "$TagLine\n";
                $count=$count+1;
        }
}

