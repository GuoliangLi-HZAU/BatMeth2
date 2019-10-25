#!/usr/bin/perl
$File=$ARGV[0];
open(FILE,$File) or die;
open(OUTPUT,">".$File.".location") or die;
print OUTPUT "0\n";
$Line=<FILE>;
while($Line=<FILE>)
{
@Des=split(/ /,$Line);
@Count=split(/ /,<FILE>);
print OUTPUT @Des[1];
print OUTPUT @Count[1]."\n";
}
