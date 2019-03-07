#use feature ':5.10';

$file = $ARGV[0];

open( INFO, $file );    # Open the file
print 0;
while ( $buff = <INFO> ) {
	if($buff =~ m/>/g){
		$buff =~ s/>//;
		print trim($count) . "\n";
		print trim($buff) . "\n";
		$count = 0;
	}else{
		$count += length(trim($buff));

	}
}
print trim($count) . "\n";
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
