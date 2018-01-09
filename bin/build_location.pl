#use feature ':5.10';
#Bowtie counter

$file = $ARGV[0];

open( INFO, $file );    # Open the file
while ( $buff = <INFO> ) {
	if($buff =~ m/>/g){
		print trim($count) . "\n";
		print trim($buff) . "\n";
		$count = 0;
	}else{
		$count += length(trim($buff));

	}


}
close(INFO);

sub trim
{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	$string =~ s/\n//g;
	$string =~ s/\r//g;
	return $string;
}
