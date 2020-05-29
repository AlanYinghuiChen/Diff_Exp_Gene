
open IN, "<test.txt";
while(<IN>){
	chomp;
	my @bam_info = split /\s+/;
	my $CIGAR = $bam_info[5];
	print "$CIGAR\n";
	@CIGAR_M = ($CIGAR =~ /(\d+)M/g);
	foreach my $num (@CIGAR_M){
		print "$num\n";
	}
}
close IN;
