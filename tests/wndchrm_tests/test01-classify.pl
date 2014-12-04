#!/usr/bin/perl -w
# Tests classification - exact sample order, exact marginal probabilities, scale factors and interpolated values.
# Also tests exact final classification accuracy.
use strict;
use warnings;
use FindBin;
use lib $FindBin::Bin;
use TestUtil;

# Forward decarations
sub parseClassifications;
sub compareClassifications;
sub getExpectedClassifications;
sub main;

# Run program
main();

# ======================================================================
sub main {

	TestUtil::exit_fail("Please specify an executable\n") unless $ARGV[0];

	my $ex = $ARGV[0];
	my $path = TestUtil::getTestPath();

	my $cmd = "$ex classify -l $path/test-l.fit $path/test-l.fit";
	print "running command:\n$cmd\n";
	my @res = `$cmd`;

	my @lines = getExpectedClassifications();
	my $expected_classifications = parseClassifications( \@lines );
	my $reported_classifications = parseClassifications( \@res );
	my $compRes;
	if (! ($compRes = compareClassifications($reported_classifications,$expected_classifications)) ) {
		print "Passed: Exact comparison (order, scale, probs., values) for ".scalar(@$expected_classifications)." classifications\n";
	} else {
		print @res;
		TestUtil::exit_fail("FAILED: $compRes\n");
	}

	my $accuracy;
	my $expected = 1.000000;
	foreach (@res) {
		$accuracy = $1 if $_ =~ /^Average accuracy.*(\d\.\d+)$/
	}
	TestUtil::exit_fail("FAILED - Accuracy not reported\n")
		unless $accuracy and TestUtil::is_numeric ($accuracy);

	if (abs ($expected - $accuracy) < TestUtil::FLT_EPSILON) {
		TestUtil::exit_pass("passed - accuracy = $accuracy, expected = $expected\n");
	} else {
		print @res;
		TestUtil::exit_fail("FAILED - accuracy = $accuracy, expected = $expected\n");
	}
}

# ======================================================================
sub parseClassifications {
	my $lines = shift;
	my @classifications;
	foreach my $line (@$lines) {
		push (@classifications, [ split (/\s+/,$line) ]) if $line =~ /[24]cell\/t\d+_s\d+_c\d+_ij.tif/;
	}
	return (\@classifications);
}

# ======================================================================
sub printClassifications {
	my $classifications = shift;

	foreach my $classification (@$classifications) {
		printf "Image [%s], scale: %g, p1: %g, p2: %g, act: [%s], pred: [%s], val: %g\n",
			$classification->[0], $classification->[1], $classification->[2], $classification->[3], 
			$classification->[4], $classification->[5], $classification->[6];
	}
}

# ======================================================================
sub compareClassifications {
	my ($reported_classifications,$expected_classifications) = (shift,shift);
	my $indx=0;

	return ("Number of reported classifications (". (1+$#{ $reported_classifications }) .") don't match number expected (". (1+$#{ $expected_classifications } ).")\n")
		unless $#{ $reported_classifications } == $#{ $expected_classifications };
	for ($indx = 0; $indx < $#{ $expected_classifications }; $indx++) {
		my $expected = $expected_classifications->[$indx];
		my $reported = $reported_classifications->[$indx];
		
		return ("Image names don't match: Expected [".$expected->[0]."], reported: [".$reported->[0]."]" )
			unless $reported->[0] eq $expected->[0];
		return ( "Scale factors don't match for image [".$reported->[0]."]: Expected (".$expected->[1]."), reported: (".$reported->[1].")" )
			unless (abs ($reported->[1] - $expected->[1]) < TestUtil::FLT_EPSILON);
		return ( "Marg. probs. don't match for image [".$reported->[0]."]: Expected (".$expected->[2].",".$expected->[3]."), reported: (".$reported->[2].",".$reported->[3].")" )
			unless (abs ($reported->[2] - $expected->[2]) < TestUtil::FLT_EPSILON && abs ($reported->[3] - $expected->[3]) < TestUtil::FLT_EPSILON);
		return ( "Interpolated values don't match for image [".$reported->[0]."]: Expected [".$expected->[6]."], reported: [".$reported->[6]."]" )
			unless (abs ($reported->[6] - $expected->[6]) < TestUtil::FLT_EPSILON);
	}
	return ("");
}

# ======================================================================
sub getExpectedClassifications  {
my $expected_out = <<END;
train/2cell/t54_s12_c04_ij.tif	1.36e-26	0.970	0.030	2cell	2cell	2.060
train/2cell/t54_s12_c06_ij.tif	1.48e-26	0.992	0.008	2cell	2cell	2.016
train/2cell/t54_s12_c08_ij.tif	4.41e-27	0.813	0.187	2cell	2cell	2.374
train/2cell/t54_s12_c11_ij.tif	8.7e-27	0.981	0.019	2cell	2cell	2.038
train/2cell/t55_s03_c03_ij.tif	8.83e-27	0.913	0.087	2cell	2cell	2.174
train/2cell/t55_s05_c09_ij.tif	5.75e-27	0.973	0.027	2cell	2cell	2.054
train/2cell/t55_s05_c12_ij.tif	7.05e-27	0.993	0.007	2cell	2cell	2.015
train/2cell/t55_s10_c11_ij.tif	1.25e-27	0.926	0.074	2cell	2cell	2.148
train/2cell/t60_s02_c10_ij.tif	5.95e-28	0.846	0.154	2cell	2cell	2.308
train/2cell/t60_s02_c11_ij.tif	2.6e-27	0.970	0.030	2cell	2cell	2.060
train/4cell/t130_s05_c06_ij.tif	3.8e-26	0.001	0.999	4cell	4cell	3.997
train/4cell/t130_s05_c09_ij.tif	2.48e-26	0.014	0.986	4cell	4cell	3.973
train/4cell/t130_s07_c08_ij.tif	1.45e-26	0.009	0.991	4cell	4cell	3.982
train/4cell/t130_s07_c11_ij.tif	1.76e-26	0.040	0.960	4cell	4cell	3.920
train/4cell/t130_s08_c05_ij.tif	4.23e-26	0.001	0.999	4cell	4cell	3.998
train/4cell/t130_s09_c07_ij.tif	2.54e-26	0.010	0.990	4cell	4cell	3.981
train/4cell/t130_s10_c01_ij.tif	3.11e-26	0.006	0.994	4cell	4cell	3.988
train/4cell/t138_s03_c02_ij.tif	1.47e-26	0.004	0.996	4cell	4cell	3.991
train/4cell/t140_s01_c10_ij.tif	2.42e-26	0.020	0.980	4cell	4cell	3.959
train/4cell/t140_s04_c01_ij.tif	4.48e-26	0.011	0.989	4cell	4cell	3.978
END

	return ( split (/\n/,$expected_out) );
}
