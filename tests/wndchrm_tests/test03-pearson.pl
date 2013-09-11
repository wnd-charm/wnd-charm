#!/usr/bin/perl -w
# Tests classification - exact sample order, exact marginal probabilities, scale factors and interpolated values.
# Also tests exact final classification accuracy.
use strict;
use warnings;
use FindBin;
use lib $FindBin::Bin;
use TestUtil;

TestUtil::exit_fail("Please specify an executable\n") unless $ARGV[0];

my $ex = $ARGV[0];
my $path = TestUtil::getTestPath();

my @res;
# First see if we can do pearson dimensionality reduction
my $do_test;
@res = `$ex`;
foreach (@res) {
	$do_test = 1 if $_ =~ /^\s+-FPearson/;
}
TestUtil::exit_pass("Passed - Executable doesn't have -FPearson option\n") unless $do_test;



@res = `$ex classify -l -FPearson $path/test-l.fit $path/test-l.fit`;



my @lines = getExpectedClassifications();
my $expected_classifications = parseClassifications (\@lines);
my $reported_classifications = parseClassifications (\@res);
my $compRes;
if (! ($compRes = compareClassifications($reported_classifications,$expected_classifications)) ) {
	print "Passed - Exact comparison (order, scale, probs., values) for ".scalar(@$expected_classifications)." classifications\n";
} else {
	print @res;
	TestUtil::exit_fail("FAILED: $compRes\n");
}


my $accuracy;
my $expected = 0.933333;
foreach (@res) {
	$accuracy = $1 if $_ =~ /^Accuracy:\s+(.+)$/;
}
TestUtil::exit_fail("FAILED - Accuracy not reported\n")
	unless $accuracy and TestUtil::is_numeric ($accuracy);

if (abs ($expected - $accuracy) < TestUtil::FLT_EPSILON) {
	TestUtil::exit_pass("Passed - accuracy = $accuracy, expected = $expected\n");
} else {
	print @res;
	TestUtil::exit_fail("FAILED - accuracy = $accuracy, expected = $expected\n");
}

sub parseClassifications {
	my $lines = shift;
	my @classifications;
	foreach my $line (@$lines) {
		push (@classifications, [ split (/\s+/,$line) ]) if $line =~ /^[24]cell\/t\d+_s\d+_c\d+_ij.tif/;
	}
	return (\@classifications);
}

sub printClassifications {
	my $classifications = shift;

	foreach my $classification (@$classifications) {
		printf "Image [%s], scale: %g, p1: %g, p2: %g, act: [%s], pred: [%s], val: %g\n",
			$classification->[0], $classification->[1], $classification->[2], $classification->[3], 
			$classification->[4], $classification->[5], $classification->[6];
	}
}

sub compareClassifications {
	my ($reported_classifications,$expected_classifications) = (shift,shift);
	my $indx=0;
	
	return ("Number of reported classifications (".scalar(@$reported_classifications).") don't match number expected (".scalar(@$expected_classifications).")\n")
		unless scalar(@$reported_classifications) == scalar(@$expected_classifications);
	for ($indx = 0; $indx < scalar(@$expected_classifications); $indx++) {
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


sub getExpectedClassifications  {
my $expected_out = <<END;
2cell/t28_s02_c10_ij.tif	9.76e-23	0.744	0.256	2cell	2cell	2.511
2cell/t29_s03_c03_ij.tif	8.67e-24	0.693	0.307	2cell	2cell	2.614
2cell/t29_s09_c02_ij.tif	1.07e-22	0.829	0.171	2cell	2cell	2.343
2cell/t32_s04_c10_ij.tif	1.77e-23	0.647	0.353	2cell	2cell	2.707
2cell/t33_s09_c01_ij.tif	6.67e-23	0.497	0.503	2cell	4cell	3.006
2cell/t33_s09_c08_ij.tif	2.83e-23	0.790	0.210	2cell	2cell	2.419
2cell/t34_s03_c12_ij.tif	1.63e-23	0.685	0.315	2cell	2cell	2.631
2cell/t35_s04_c01_ij.tif	2.5e-23	0.822	0.178	2cell	2cell	2.356
2cell/t36_s03_c11_ij.tif	5.14e-23	0.713	0.287	2cell	2cell	2.575
2cell/t36_s10_c08_ij.tif	3.96e-23	0.873	0.127	2cell	2cell	2.254
2cell/t38_s02_c04_ij.tif	2.44e-23	0.720	0.280	2cell	2cell	2.561
2cell/t40_s04_c04_ij.tif	3.7e-23	0.416	0.584	2cell	4cell	3.167
2cell/t40_s07_c08_ij.tif	4.9e-23	0.568	0.432	2cell	2cell	2.864
2cell/t40_s10_c12_ij.tif	9.5e-23	0.717	0.283	2cell	2cell	2.566
2cell/t40_s9_c09_ij.tif	2.23e-23	0.728	0.272	2cell	2cell	2.544
2cell/t43_s03_c03_ij.tif	2.55e-23	0.740	0.260	2cell	2cell	2.520
2cell/t43_s03_c12_ij.tif	5.47e-23	0.562	0.438	2cell	2cell	2.876
2cell/t45_s01_c10_ij.tif	9.72e-23	0.587	0.413	2cell	2cell	2.826
2cell/t45_s06_c03_ij.tif	4.99e-23	0.911	0.089	2cell	2cell	2.178
2cell/t45_s06_c06_ij.tif	1.48e-23	0.369	0.631	2cell	4cell	3.262
2cell/t45_s06_c12_ij.tif	7.64e-23	0.704	0.296	2cell	2cell	2.591
2cell/t45_s07_c05_ij.tif	7.32e-23	0.808	0.192	2cell	2cell	2.384
2cell/t45_s07_c11_ij.tif	4.04e-23	0.909	0.091	2cell	2cell	2.182
2cell/t45_s11_c10_ij.tif	4.62e-23	0.424	0.576	2cell	4cell	3.151
2cell/t45_s13_c04_ij.tif	1.37e-23	0.589	0.411	2cell	2cell	2.822
2cell/t47_s08_c07_ij.tif	1.41e-23	0.785	0.215	2cell	2cell	2.430
2cell/t47_s08_c08_ij.tif	1.8e-23	0.789	0.211	2cell	2cell	2.421
2cell/t47_s8_c04_ij.tif	1.05e-22	0.755	0.245	2cell	2cell	2.490
2cell/t48_s11_c03_ij.tif	2.94e-23	0.711	0.289	2cell	2cell	2.577
2cell/t50_s04_c05_ij.tif	2.43e-23	0.782	0.218	2cell	2cell	2.436
2cell/t50_s05_c11_ij.tif	3.15e-23	0.737	0.263	2cell	2cell	2.526
2cell/t50_s10_c05_ij.tif	5.68e-23	0.842	0.158	2cell	2cell	2.316
2cell/t52_s3_c10_ij.tif	2.81e-23	0.717	0.283	2cell	2cell	2.566
2cell/t52_s9_c10_ij.tif	4.45e-23	0.818	0.182	2cell	2cell	2.364
2cell/t54_s12_c03_ij.tif	4.04e-23	0.819	0.181	2cell	2cell	2.361
2cell/t54_s12_c04_ij.tif	9.14e-24	0.637	0.363	2cell	2cell	2.727
2cell/t54_s12_c06_ij.tif	1.21e-23	0.849	0.151	2cell	2cell	2.302
2cell/t54_s12_c08_ij.tif	3.64e-23	0.877	0.123	2cell	2cell	2.245
2cell/t54_s12_c11_ij.tif	2.83e-23	0.813	0.187	2cell	2cell	2.374
2cell/t55_s03_c03_ij.tif	7.67e-23	0.705	0.295	2cell	2cell	2.591
2cell/t55_s05_c09_ij.tif	6.71e-23	0.479	0.521	2cell	4cell	3.042
2cell/t55_s05_c12_ij.tif	1.89e-23	0.446	0.554	2cell	4cell	3.107
2cell/t55_s10_c11_ij.tif	2.38e-23	0.766	0.234	2cell	2cell	2.468
2cell/t60_s02_c10_ij.tif	2.03e-24	0.527	0.473	2cell	2cell	2.946
2cell/t60_s02_c11_ij.tif	2.57e-23	0.703	0.297	2cell	2cell	2.595
4cell/t103_s12_c11_ij.tif	7.27e-23	0.223	0.777	4cell	4cell	3.554
4cell/t106_s01_c09_ij.tif	5.84e-23	0.264	0.736	4cell	4cell	3.472
4cell/t108_s03_c10_ij.tif	9.36e-23	0.121	0.879	4cell	4cell	3.759
4cell/t108_s06_c04_ij.tif	3.42e-23	0.166	0.834	4cell	4cell	3.668
4cell/t108_s06_c09_ij.tif	7.78e-23	0.162	0.838	4cell	4cell	3.677
4cell/t109_s11_c09_ij.tif	5.4e-23	0.073	0.927	4cell	4cell	3.855
4cell/t110_s03_c12_ij.tif	3.33e-23	0.409	0.591	4cell	4cell	3.182
4cell/t110_s10_c12_ij.tif	9.99e-23	0.132	0.868	4cell	4cell	3.736
4cell/t113_s09_c10_ij.tif	1.01e-22	0.155	0.845	4cell	4cell	3.689
4cell/t114_s01_c06_ij.tif	5.21e-23	0.277	0.723	4cell	4cell	3.445
4cell/t114_s08_c07_ij.tif	8.25e-23	0.189	0.811	4cell	4cell	3.623
4cell/t115_s02_c12_ij.tif	7.13e-23	0.115	0.885	4cell	4cell	3.769
4cell/t117_s05_c05_ij.tif	1.38e-22	0.199	0.801	4cell	4cell	3.601
4cell/t118_s04_c12_ij.tif	8.96e-23	0.110	0.890	4cell	4cell	3.780
4cell/t120_s06_c05_ij.tif	1.41e-22	0.149	0.851	4cell	4cell	3.703
4cell/t120_s06_c06_ij.tif	3.28e-23	0.092	0.908	4cell	4cell	3.817
4cell/t120_s06_c12_ij.tif	1e-22	0.214	0.786	4cell	4cell	3.573
4cell/t120_s12_c02_ij.tif	4e-23	0.224	0.776	4cell	4cell	3.553
4cell/t120_s12_c06_ij.tif	1.02e-22	0.072	0.928	4cell	4cell	3.856
4cell/t120_s12_c07_ij.tif	1.18e-22	0.090	0.910	4cell	4cell	3.819
4cell/t123_s04_c02_ij.tif	3.85e-23	0.202	0.798	4cell	4cell	3.597
4cell/t128_s03_c11_ij.tif	9.33e-23	0.122	0.878	4cell	4cell	3.756
4cell/t130_s02_c01_ij.tif	1.52e-22	0.081	0.919	4cell	4cell	3.838
4cell/t130_s05_c06_ij.tif	1.14e-22	0.153	0.847	4cell	4cell	3.694
4cell/t130_s05_c09_ij.tif	7.64e-23	0.072	0.928	4cell	4cell	3.856
4cell/t130_s07_c08_ij.tif	3.69e-23	0.284	0.716	4cell	4cell	3.432
4cell/t130_s07_c11_ij.tif	8.18e-23	0.261	0.739	4cell	4cell	3.478
4cell/t130_s08_c05_ij.tif	1.35e-22	0.085	0.915	4cell	4cell	3.830
4cell/t130_s09_c07_ij.tif	5.32e-23	0.158	0.842	4cell	4cell	3.684
4cell/t130_s10_c01_ij.tif	9.23e-23	0.264	0.736	4cell	4cell	3.472
4cell/t138_s03_c02_ij.tif	7.05e-23	0.094	0.906	4cell	4cell	3.812
4cell/t140_s01_c10_ij.tif	4.66e-23	0.130	0.870	4cell	4cell	3.741
4cell/t140_s04_c01_ij.tif	1.31e-22	0.096	0.904	4cell	4cell	3.808
4cell/t1_s01_c05_ij.tif	6.97e-23	0.069	0.931	4cell	4cell	3.861
4cell/t1_s01_c10_ij.tif	6.42e-23	0.084	0.916	4cell	4cell	3.832
4cell/t1_s02_c11_ij.tif	7.53e-23	0.256	0.744	4cell	4cell	3.487
4cell/t1_s03_c08_ij.tif	2.44e-23	0.420	0.580	4cell	4cell	3.159
4cell/t1_s04_c04_ij.tif	8.97e-23	0.117	0.883	4cell	4cell	3.767
4cell/t1_s05_c06_ij.tif	3.43e-23	0.231	0.769	4cell	4cell	3.539
4cell/t1_s07_c08_ij.tif	7.17e-23	0.175	0.825	4cell	4cell	3.650
4cell/t1_s09_c01_ij.tif	4.43e-23	0.135	0.865	4cell	4cell	3.729
4cell/t1_s09_c11_ij.tif	7.13e-23	0.219	0.781	4cell	4cell	3.562
4cell/t1_s10_c08_ij.tif	1.08e-22	0.143	0.857	4cell	4cell	3.713
4cell/t1_s13_c04_ij.tif	9.81e-23	0.144	0.856	4cell	4cell	3.712
4cell/t1_s14_c06_ij.tif	1.81e-22	0.065	0.935	4cell	4cell	3.870
END

	return ( split (/\n/,$expected_out) );
}
