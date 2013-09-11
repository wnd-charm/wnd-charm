#!/usr/bin/perl -w
# Tests classification - exact sample order, exact marginal probabilities, scale factors and interpolated values.
# Also tests exact final classification accuracy.
use strict;
use warnings;
use FindBin;
use lib $FindBin::Bin;
use TestUtil;
use File::Find;

TestUtil::exit_fail("Please specify an executable\n") unless $ARGV[0];

my $ex = $ARGV[0];
my $path = TestUtil::getTestPath();

print "removing sigs in $path/images...\n";
find(\&do_entry, "$path/images");

print "calculating sigs...\n";
`cd $path; $ex train -l images test-images-sigs-l.fit`;

print "comparing to unbalanced.fit\n";
my @res = `diff $path/test-images-sigs-l.fit $path/unbalanced.fit`;
chomp foreach (@res);

if (length (@res) < 2) {
	print "Passed - No differences found\n";
	unlink ("$path/test-images-sigs-l.fit");
} else {
	print "Failed - differences found in sigs:\ndiff $path/test-images-sigs-l.fit $path/unbalanced.fit\n";
}

sub do_entry {
	/\.sig$/ and unlink;
}
