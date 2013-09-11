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
train/2cell/t54_s12_c04_ij.tif	9.3e-27	0.966	0.034	2cell	2cell	2.068
train/2cell/t54_s12_c06_ij.tif	8.7e-27	0.989	0.011	2cell	2cell	2.021
train/2cell/t54_s12_c08_ij.tif	4.54e-27	0.838	0.162	2cell	2cell	2.324
train/2cell/t54_s12_c11_ij.tif	9.88e-27	0.983	0.017	2cell	2cell	2.033
train/2cell/t55_s03_c03_ij.tif	1.01e-26	0.922	0.078	2cell	2cell	2.156
train/2cell/t55_s05_c09_ij.tif	6.31e-27	0.976	0.024	2cell	2cell	2.048
train/2cell/t55_s05_c12_ij.tif	6.04e-27	0.993	0.007	2cell	2cell	2.014
train/2cell/t55_s10_c11_ij.tif	1.15e-27	0.912	0.088	2cell	2cell	2.176
train/2cell/t60_s02_c10_ij.tif	6.64e-28	0.878	0.122	2cell	2cell	2.244
train/2cell/t60_s02_c11_ij.tif	2.18e-27	0.969	0.031	2cell	2cell	2.063
train/4cell/t130_s05_c06_ij.tif	2.39e-26	0.002	0.998	4cell	4cell	3.996
train/4cell/t130_s05_c09_ij.tif	2.29e-26	0.016	0.984	4cell	4cell	3.967
train/4cell/t130_s07_c08_ij.tif	1.15e-26	0.009	0.991	4cell	4cell	3.983
train/4cell/t130_s07_c11_ij.tif	1.49e-26	0.032	0.968	4cell	4cell	3.936
train/4cell/t130_s08_c05_ij.tif	2.5e-26	0.001	0.999	4cell	4cell	3.997
train/4cell/t130_s09_c07_ij.tif	3.32e-26	0.007	0.993	4cell	4cell	3.986
train/4cell/t130_s10_c01_ij.tif	2.74e-26	0.008	0.992	4cell	4cell	3.983
train/4cell/t138_s03_c02_ij.tif	1.35e-26	0.006	0.994	4cell	4cell	3.988
train/4cell/t140_s01_c10_ij.tif	1.78e-26	0.022	0.978	4cell	4cell	3.956
train/4cell/t140_s04_c01_ij.tif	3.74e-26	0.015	0.985	4cell	4cell	3.970
END

	return ( split (/\n/,$expected_out) );
}
