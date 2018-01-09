use strict;
die "perl $0 snp.txt methylation.txt methylation_filtered.txt methylation_SNP.txt\n" unless(@ARGV==4);
my $snp=shift;
my $cpg=shift;
my $cgl=shift;
my $cgf=shift;
my %hash;

open SNP, $snp or die $!;
open CG,$cpg or die $!;
open CGL, ">$cgl" or die $!;
open CGF, ">$cgf" or die $!;

while(<SNP>){
	chomp;
	my @a=split;
	if($a[3] eq "C" && $a[4] eq "T"){
		$hash{$a[0]}{$a[1]}=1;
	}	
	if($a[3] eq "G" && $a[4] eq "A"){
		my $pos=$a[1]-1;
		$hash{$a[0]}{$pos}=1;

	}

}
close SNP;

while(<CG>){
	chomp;
	my @a=split;
	if(exists($hash{$a[0]}{$a[1]})){
		print CGF $_."\n";
	}else{
		print CGL $_."\n";
	}
	
}
close CG;
