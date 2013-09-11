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



my @lines = getExpectedClassifications();
my ($expected_2as4,$expected_4as2,$expected_accuracies,$expected_avg) = parseClassifications (\@lines);
my $cmd = "$ex test -ls -i25 -j20 -n100 $path/test-l.fit";

# Note that we only need to pass this test once
my $compRes;
my @res;
my $tries = 2;
while ($tries) {
	@res = `$cmd`;
	my ($reported_2as4,$reported_4as2,$reported_accuracies,$reported_avg) = parseClassifications (\@res);

	print @res and TestUtil::exit_fail("FAILED: Number of reported classifications does not match expected:\n".
		"    2 as 4 expected: ".scalar (@$expected_2as4).", got ".scalar (@$reported_2as4)."\n".
		"    4 as 2 expected: ".scalar (@$expected_4as2).", got ".scalar (@$reported_4as2)."\n".
		"accuracies expected: ".scalar (@$expected_accuracies).", got ".scalar (@$reported_accuracies)."\n")
		unless
			scalar (@$reported_2as4) == scalar (@$expected_2as4) and
			scalar (@$reported_4as2) == scalar (@$expected_4as2) and
			scalar (@$reported_accuracies) == scalar (@$expected_accuracies);

	if (! ($compRes = compareClassifications($reported_2as4,$reported_4as2,$reported_accuracies,$reported_avg,$expected_2as4,$expected_4as2,$expected_accuracies,$expected_avg)) ) {
		print "Passed - Distribution comparison (2 as 4, 4 as 2, accuracies) for ".scalar(@$reported_accuracies)." classifications.\n";
		last;
	}
	$tries--;
}
if ($compRes) {
	print @res;
	TestUtil::exit_fail("FAILED: $compRes\n");
}


sub parseClassifications {
	my $lines = shift;
	my @class2as4;
	my @class4as2;
	my @accuracies;
	my $avg_accuracy;
	foreach my $line (@$lines) {
		push (@class2as4, $1) if $line =~ /^\s+2cell\s+\d\.\d+\s+(\d\.\d+)/;
		push (@class4as2, $1) if $line =~ /^\s+4cell\s+(\d\.\d+)\s+\d\.\d+/;
		push (@accuracies, $1) if $line =~ /^Accuracy:\s+(\d\.\d+)\s*$/;
		$avg_accuracy = $1 if $line =~ /^Average accuracy.+(\d\.\d+)$/;
	}
	
	return (\@class2as4,\@class4as2,\@accuracies,$avg_accuracy);
}

sub printClassifications {
	my ($class_2as4,$class_4as2,$class_accuracies,$reported_avg) = (shift,shift,shift,shift);
	my $i;
	
	for ($i=0; $i < scalar (@$class_2as4); $i++) {
		printf "%.4f\t%.4f\t%.4f\n",$class_2as4->[$i],$class_4as2->[$i],$class_accuracies->[$i];
	}
	printf "Average: %g\n",$reported_avg;
}

sub compareClassifications {
	my ($reported_2as4,$reported_4as2,$reported_accuracies,$reported_avg,$expected_2as4,$expected_4as2,$expected_accuracies,$expected_avg) = 
		(shift,shift,shift,shift,shift,shift,shift,shift);

	my ($exp_acc_mean,$exp_acc_stddev);
	my ($rep_acc_mean,$rep_acc_stddev);
	my ($t_prob,$f_prob);

 	$t_prob = TestUtil::t_test($expected_2as4,$reported_2as4);
 	return "Class 2 as 4 distributions differ (P=$t_prob)" unless $t_prob > 0.15;
	$f_prob = TestUtil::f_test($expected_2as4,$reported_2as4,99.9);
	return "Class 2 as 4 distributions have different variances (P<99.9%)" unless $f_prob;

	$t_prob = TestUtil::t_test($expected_4as2,$reported_4as2);
	return "Class 4 as 2 distributions differ (P=$t_prob)" unless $t_prob > 0.15;
	$f_prob = TestUtil::f_test($expected_4as2,$reported_4as2,99.9);
	return "Class 4 as 2 distributions have different variances (P<99.9%)" unless $f_prob;

	$t_prob = TestUtil::t_test($expected_accuracies,$reported_accuracies);
	return "Accuracy distributions differ (P=$t_prob)" unless $t_prob > 0.15;
	$f_prob = TestUtil::f_test($expected_accuracies,$reported_accuracies,99.9);
	return "Accuracy distributions have different variances (P<99.9%)" unless $f_prob;
	
	my $mean = TestUtil::getMean($reported_accuracies);
	$reported_avg = sprintf "%.5f",$reported_avg;
	$mean = sprintf "%.5f",$mean;
	return "Average accuracy ($mean) does not match reported average ($reported_avg)" unless abs ($mean-$reported_avg) < TestUtil::FLT_EPSILON;

	return ("");
}


sub getExpectedClassifications  {
my $expected_out = <<END;
           2cell         1.00000         0.20901
           4cell         0.13390         1.00000
Accuracy: 0.925000 
           2cell         1.00000         0.24596
           4cell         0.21109         1.00000
Accuracy: 0.900000 
           2cell         1.00000         0.42783
           4cell         0.15046         1.00000
Accuracy: 0.825000 
           2cell         1.00000         0.36165
           4cell         0.15460         1.00000
Accuracy: 0.900000 
           2cell         1.00000         0.31977
           4cell         0.13163         1.00000
Accuracy: 0.925000 
           2cell         1.00000         0.26521
           4cell         0.08409         1.00000
Accuracy: 0.950000 
           2cell         1.00000         0.58745
           4cell         0.07448         1.00000
Accuracy: 0.825000 
           2cell         1.00000         0.19958
           4cell         0.22819         1.00000
Accuracy: 0.900000 
           2cell         1.00000         0.33258
           4cell         0.12658         1.00000
Accuracy: 0.925000 
           2cell         1.00000         0.26141
           4cell         0.15074         1.00000
Accuracy: 0.950000 
           2cell         1.00000         0.29062
           4cell         0.17098         1.00000
Accuracy: 0.900000 
           2cell         1.00000         0.35343
           4cell         0.13038         1.00000
Accuracy: 0.900000 
           2cell         1.00000         0.21251
           4cell         0.23752         1.00000
Accuracy: 0.950000 
           2cell         1.00000         0.38034
           4cell         0.20834         1.00000
Accuracy: 0.825000 
           2cell         1.00000         0.27721
           4cell         0.19663         1.00000
Accuracy: 0.925000 
           2cell         1.00000         0.23689
           4cell         0.21488         1.00000
Accuracy: 0.925000 
           2cell         1.00000         0.38862
           4cell         0.13017         1.00000
Accuracy: 0.900000 
           2cell         1.00000         0.36789
           4cell         0.19717         1.00000
Accuracy: 0.925000 
           2cell         1.00000         0.32667
           4cell         0.18473         1.00000
Accuracy: 0.900000 
           2cell         1.00000         0.26426
           4cell         0.14430         1.00000
Accuracy: 0.950000 
           2cell         1.00000         0.18515
           4cell         0.15817         1.00000
Accuracy: 0.925000 
           2cell         1.00000         0.27999
           4cell         0.18147         1.00000
Accuracy: 0.950000 
           2cell         1.00000         0.39906
           4cell         0.13770         1.00000
Accuracy: 0.875000 
           2cell         1.00000         0.26610
           4cell         0.19597         1.00000
Accuracy: 0.900000 
           2cell         1.00000         0.47761
           4cell         0.14858         1.00000
Accuracy: 0.875000 
           2cell         1.00000         0.29186
           4cell         0.18429         1.00000
Accuracy: 0.950000 
           2cell         1.00000         0.24801
           4cell         0.16100         1.00000
Accuracy: 0.975000 
           2cell         1.00000         0.28238
           4cell         0.13741         1.00000
Accuracy: 0.925000 
           2cell         1.00000         0.44728
           4cell         0.13797         1.00000
Accuracy: 0.925000 
           2cell         1.00000         0.31669
           4cell         0.20222         1.00000
Accuracy: 0.900000 
           2cell         1.00000         0.17476
           4cell         0.19972         1.00000
Accuracy: 0.925000 
           2cell         1.00000         0.27767
           4cell         0.15635         1.00000
Accuracy: 0.950000 
           2cell         1.00000         0.35524
           4cell         0.24964         1.00000
Accuracy: 0.900000 
           2cell         1.00000         0.47075
           4cell         0.07854         1.00000
Accuracy: 0.875000 
           2cell         1.00000         0.47303
           4cell         0.10439         1.00000
Accuracy: 0.900000 
           2cell         1.00000         0.24572
           4cell         0.21013         1.00000
Accuracy: 0.900000 
           2cell         1.00000         0.50568
           4cell         0.17726         1.00000
Accuracy: 0.800000 
           2cell         1.00000         0.18682
           4cell         0.25215         1.00000
Accuracy: 0.900000 
           2cell         1.00000         0.21320
           4cell         0.20622         1.00000
Accuracy: 0.950000 
           2cell         1.00000         0.46940
           4cell         0.14500         1.00000
Accuracy: 0.900000 
           2cell         1.00000         0.26858
           4cell         0.15531         1.00000
Accuracy: 0.900000 
           2cell         1.00000         0.28011
           4cell         0.16511         1.00000
Accuracy: 0.950000 
           2cell         1.00000         0.41583
           4cell         0.08533         1.00000
Accuracy: 0.875000 
           2cell         1.00000         0.23566
           4cell         0.17224         1.00000
Accuracy: 0.950000 
           2cell         1.00000         0.41142
           4cell         0.11180         1.00000
Accuracy: 0.925000 
           2cell         1.00000         0.25684
           4cell         0.15449         1.00000
Accuracy: 0.925000 
           2cell         1.00000         0.30701
           4cell         0.12766         1.00000
Accuracy: 0.925000 
           2cell         1.00000         0.42175
           4cell         0.26525         1.00000
Accuracy: 0.900000 
           2cell         1.00000         0.27203
           4cell         0.19141         1.00000
Accuracy: 0.925000 
           2cell         1.00000         0.30840
           4cell         0.18071         1.00000
Accuracy: 0.850000 
           2cell         1.00000         0.24484
           4cell         0.18728         1.00000
Accuracy: 0.950000 
           2cell         1.00000         0.33535
           4cell         0.13141         1.00000
Accuracy: 0.925000 
           2cell         1.00000         0.34816
           4cell         0.10266         1.00000
Accuracy: 0.925000 
           2cell         1.00000         0.20442
           4cell         0.20175         1.00000
Accuracy: 0.925000 
           2cell         1.00000         0.27654
           4cell         0.16989         1.00000
Accuracy: 0.925000 
           2cell         1.00000         0.44126
           4cell         0.18899         1.00000
Accuracy: 0.850000 
           2cell         1.00000         0.40318
           4cell         0.09534         1.00000
Accuracy: 0.925000 
           2cell         1.00000         0.41702
           4cell         0.15116         1.00000
Accuracy: 0.875000 
           2cell         1.00000         0.29559
           4cell         0.23654         1.00000
Accuracy: 0.900000 
           2cell         1.00000         0.43399
           4cell         0.42457         1.00000
Accuracy: 0.725000 
           2cell         1.00000         0.19734
           4cell         0.17884         1.00000
Accuracy: 0.975000 
           2cell         1.00000         0.31793
           4cell         0.12480         1.00000
Accuracy: 0.950000 
           2cell         1.00000         0.18063
           4cell         0.14159         1.00000
Accuracy: 0.950000 
           2cell         1.00000         0.44579
           4cell         0.20469         1.00000
Accuracy: 0.875000 
           2cell         1.00000         0.28509
           4cell         0.14058         1.00000
Accuracy: 0.900000 
           2cell         1.00000         0.29386
           4cell         0.17542         1.00000
Accuracy: 0.925000 
           2cell         1.00000         0.28147
           4cell         0.12523         1.00000
Accuracy: 0.950000 
           2cell         1.00000         0.34292
           4cell         0.15843         1.00000
Accuracy: 0.825000 
           2cell         1.00000         0.15929
           4cell         0.15232         1.00000
Accuracy: 0.975000 
           2cell         1.00000         0.50624
           4cell         0.15051         1.00000
Accuracy: 0.800000 
           2cell         1.00000         0.38189
           4cell         0.19792         1.00000
Accuracy: 0.825000 
           2cell         1.00000         0.27214
           4cell         0.14108         1.00000
Accuracy: 0.950000 
           2cell         1.00000         0.39863
           4cell         0.12725         1.00000
Accuracy: 0.925000 
           2cell         1.00000         0.39819
           4cell         0.22596         1.00000
Accuracy: 0.850000 
           2cell         1.00000         0.20649
           4cell         0.22765         1.00000
Accuracy: 0.900000 
           2cell         1.00000         0.24656
           4cell         0.16402         1.00000
Accuracy: 0.950000 
           2cell         1.00000         0.21365
           4cell         0.12224         1.00000
Accuracy: 0.950000 
           2cell         1.00000         0.20101
           4cell         0.16678         1.00000
Accuracy: 0.975000 
           2cell         1.00000         0.30194
           4cell         0.13499         1.00000
Accuracy: 0.975000 
           2cell         1.00000         0.39384
           4cell         0.13120         1.00000
Accuracy: 0.925000 
           2cell         1.00000         0.29752
           4cell         0.12185         1.00000
Accuracy: 0.900000 
           2cell         1.00000         0.36017
           4cell         0.29466         1.00000
Accuracy: 0.850000 
           2cell         1.00000         0.22484
           4cell         0.24882         1.00000
Accuracy: 0.875000 
           2cell         1.00000         0.30084
           4cell         0.09697         1.00000
Accuracy: 0.925000 
           2cell         1.00000         0.36194
           4cell         0.14156         1.00000
Accuracy: 0.975000 
           2cell         1.00000         0.37363
           4cell         0.16087         1.00000
Accuracy: 0.850000 
           2cell         1.00000         0.30765
           4cell         0.13928         1.00000
Accuracy: 0.875000 
           2cell         1.00000         0.44616
           4cell         0.09451         1.00000
Accuracy: 0.875000 
           2cell         1.00000         0.14636
           4cell         0.18980         1.00000
Accuracy: 0.950000 
           2cell         1.00000         0.34974
           4cell         0.22151         1.00000
Accuracy: 0.975000 
           2cell         1.00000         0.37108
           4cell         0.20578         1.00000
Accuracy: 0.850000 
           2cell         1.00000         0.12041
           4cell         0.30323         1.00000
Accuracy: 0.925000 
           2cell         1.00000         0.41630
           4cell         0.15515         1.00000
Accuracy: 0.875000 
           2cell         1.00000         0.50440
           4cell         0.16556         1.00000
Accuracy: 0.800000 
           2cell         1.00000         0.26369
           4cell         0.20520         1.00000
Accuracy: 0.900000 
           2cell         1.00000         0.17811
           4cell         0.19451         1.00000
Accuracy: 0.975000 
           2cell         1.00000         0.17973
           4cell         0.15260         1.00000
Accuracy: 1.000000 
           2cell         1.00000         0.35401
           4cell         0.16834         1.00000
Accuracy: 1.000000 
           2cell         1.00000         0.19272
           4cell         0.11633         1.00000
Accuracy: 0.950000 
           2cell         1.00000         0.33216
           4cell         0.08070         1.00000
Accuracy: 0.900000 
Average accuracy (100 splits): 0.909500
END

	return ( split (/\n/,$expected_out) );
}
