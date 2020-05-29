use strict;
use FindBin '$Bin';
use Getopt::Long;
use File::Basename;
my ($inBam,$outBam,$samtools,$Rscript,$help);

GetOptions(
	"inBam|i:s" => \$inBam,
	"outBam|o:s" => \$outBam,
	"samtools|s:s" =>\$samtools,
	"Rscript|r:s" =>\$Rscript,
	"help|?" => \$help
);

$samtools ||= "/BioII/lulab_b/chenyinghui/software/conda3/bin/samtools";
$Rscript  ||= "/BioII/lulab_b/chenyinghui/software/conda3/bin/Rscript";

my $outSam = "$outBam".".sam";
my $outBai = "$outBam".".bai";

open IN, "$samtools view -h $inBam |" or die "can't open $inBam: $!\n";
open(TEMP,">$outSam");
while(<IN>){
	if(/jM:B:c,-1/){
		next;
	}else{
		print TEMP $_;
	}
}
close IN;
close TEMP;

system("$samtools view -hb -o $outBam $outSam");
system("rm $outSam");
system("$samtools index -b $outBam > $outBai");
