#!/usr/bin/perl -w
use strict;
use warnings;
use FindBin;
package TestUtil;
use Statistics::TTest;

use constant FLT_EPSILON => 1.19209290e-07;
use constant EXIT_SUCCESS => 0;
use constant EXIT_FAILURE => -1;

sub getnum {
	use POSIX qw(strtod);
	my $str = shift;
	$str =~ s/^\s+//;
	$str =~ s/\s+$//;
	$! = 0;
	my($num, $unparsed) = strtod($str);
	if (($str eq '') || ($unparsed != 0) || $!) {
		return;
	} else {
		return $num;
	}
}
sub is_numeric { defined scalar &getnum }

sub getTestPath { return  ($FindBin::Bin); }

sub exit_pass { print shift; exit (EXIT_SUCCESS); }
sub exit_fail { print shift; exit (EXIT_FAILURE); }

sub t_test {
	my ($r1,$r2,$sig) = (shift,shift,shift);
	$sig = 99.9 unless $sig;
	my $ttest = new Statistics::TTest;
	$ttest->load_data($r1,$r2);
	$ttest->set_significance($sig);
	return $ttest->{t_prob};
}

sub f_test {
	my ($r1,$r2,$sig) = (shift,shift,shift);
	my $do_cutoff = 1 if $sig;
	$sig = 99.9 unless $sig;
	my $ttest = new Statistics::TTest;
	$ttest->load_data($r1,$r2);
	$ttest->set_significance($sig);
	return $ttest->{f_statistic} unless $do_cutoff;
	return ($ttest->{f_statistic} < $ttest->{f_cutoff});
}

sub getMean {
	my $stat = Statistics::Descriptive::Sparse->new();
	$stat->add_data(shift);
	return ($stat->mean());
}


return 1;
