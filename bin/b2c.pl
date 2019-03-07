#!/usr/bin/perl 
print "Base space 2 Color space...\n";

$FileName=$ARGV[0];
unless(open(FILE,$FileName))
{
	print "File does not Exist \n";
	exit;
}

%Conversion_Table=();
$Conversion_Table{'AA'}='A';
$Conversion_Table{'CC'}='A';
$Conversion_Table{'GG'}='A';
$Conversion_Table{'TT'}='A';

$Conversion_Table{'AC'}='C';
$Conversion_Table{'CA'}='C';
$Conversion_Table{'GT'}='C';
$Conversion_Table{'TG'}='C';

$Conversion_Table{'AG'}='G';
$Conversion_Table{'GA'}='G';
$Conversion_Table{'CT'}='G';
$Conversion_Table{'TC'}='G';

$Conversion_Table{'AT'}='T';
$Conversion_Table{'TA'}='T';
$Conversion_Table{'GC'}='T';
$Conversion_Table{'CG'}='T';

if($ARGV[1]) {open (TAGFILE,'>'.$ARGV[1]) or die;}
else {open (TAGFILE,'>'.$FileName.".cfq") or die;}
$Line_Count = 0;
$count = 0;
$Last='';
$New_Tag='';#add first char to be .
$Pass=0;
while ($TagLine = <FILE> )
{
	if (substr($TagLine,0,1) eq '>')
	{
		if (!$Pass) 
		{
			print TAGFILE $TagLine;
			$Last=getc FILE;
			print TAGFILE $Last;
			$Pass=1;
		}
		else
		{
			print TAGFILE $New_Tag."\n";
			$New_Tag='';
			print TAGFILE $TagLine;
		}
	}
	else
	{
		$TagLine=$Last.$TagLine;
		chomp $TagLine;
		$TagLine=~ tr/nNacgt/AAACGT/;
		$count=$count+1;
		#$New_Tag='';
		$Len= (length $TagLine);
		for ($i=0;$i< $Len-1 ; $i++)
		{
			$Current=substr($TagLine,$i,2);
			$New_Tag=$New_Tag.$Conversion_Table{$Current};
			if((length $New_Tag)==50)
			{
				print TAGFILE $New_Tag."\n";
				$New_Tag='';
			}
		}
		#$New_Tag=$New_Tag.$Conversion_Table{substr($TagLine,$_,2)} for 0 .. length $TagLine;
		$Last=substr($TagLine,(length $TagLine)-1,1);
	}
}
print TAGFILE $New_Tag."\n";
print $count." lines Coded... \n";
close FILE;
