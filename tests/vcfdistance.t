use strict;
use warnings;
use Test::More;
use Local::vcflib::Test;

my @vcf = split /\n/, <<'';
##fileformat=VCFv4.0
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
refseq	502	.	G	A	38553	PASS	
refseq	552	.	G	A	24044	PASS	
refseq	660	.	G	A	38553	PASS	
refseq	678	.	G	A	24044	PASS	
refseq	684	.	G	A	24044	PASS	

sub variants {
    join "\n", @vcf[0 .. $_[0] + 1]
}

my ($output, $header) = ('', <<'');
##fileformat=VCFv4.0
##INFO=<ID=BasesToClosestVariant,Number=1,Type=Integer,Description="Number of bases to the closest variant in the file.">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO

($output) = run_ok(["vcfdistance"], variants(5));
is $output, $header . <<'', "distances for 5 variant lines";
refseq	502	.	G	A	38553	PASS	BasesToClosestVariant=50;
refseq	552	.	G	A	24044	PASS	BasesToClosestVariant=50;
refseq	660	.	G	A	38553	PASS	BasesToClosestVariant=18;
refseq	678	.	G	A	24044	PASS	BasesToClosestVariant=6;
refseq	684	.	G	A	24044	PASS	BasesToClosestVariant=6;

($output) = run_ok(["vcfdistance"], variants(4));
is $output, $header . <<'', "distances for 4 variant lines";
refseq	502	.	G	A	38553	PASS	BasesToClosestVariant=50;
refseq	552	.	G	A	24044	PASS	BasesToClosestVariant=50;
refseq	660	.	G	A	38553	PASS	BasesToClosestVariant=18;
refseq	678	.	G	A	24044	PASS	BasesToClosestVariant=18;

($output) = run_ok(["vcfdistance"], variants(3));
is $output, $header . <<'', "distances for 3 variant lines";
refseq	502	.	G	A	38553	PASS	BasesToClosestVariant=50;
refseq	552	.	G	A	24044	PASS	BasesToClosestVariant=50;
refseq	660	.	G	A	38553	PASS	BasesToClosestVariant=108;

($output) = run_ok(["vcfdistance"], variants(2));
is $output, $header . <<'', "distances for 2 variant lines";
refseq	502	.	G	A	38553	PASS	BasesToClosestVariant=50;
refseq	552	.	G	A	24044	PASS	BasesToClosestVariant=50;

($output) = run_ok(["vcfdistance"], variants(1));
is $output, $header . <<'', "distances for 1 variant line";
refseq	502	.	G	A	38553	PASS	

($output) = run_ok(["vcfdistance"], variants(0));
is $output, $header, "distances for 0 variant lines";

done_testing;
