#!/usr/bin/perl
##
##
###############################################################################

use strict;
use POSIX;

if ($#ARGV < 0) {
	print STDERR "usage: $0 <scoretag> <filelist>\n";
	exit -1;
}

my @RMSBINS=(
	0.0, 0.5, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.5, 4.0, 5.0, 6.0, 9999.0
); ## 9999 is a placeholder
my $nrmsbins = $#RMSBINS;

my $TEMP = 0.1;    # scaled by 5-95% spread

# params (flag?)
my $MIN_BIN_SIZE = 4;
my $verbose = 1;
my $cap=1.0;

my $scoretag = shift @ARGV;
my @silentfiles;
my $tag = shift (@ARGV);
open (FILELIST, $tag) || die "Unable to open $tag.";
my @silentfiles = <FILELIST>;
chomp (@silentfiles);
close (FILELIST);

my $nsilent = scalar( @silentfiles );

my $SCORECAP = 0;

my $avg_descrim = 0;

my $sumvalidbins = 0;
my $sumvalidstructs = 0;

my @all_counts;
my @all_energies;
my @all_rms;

foreach my $j (0..$nsilent-1) {
	open (SIL, $silentfiles[$j]) || die "Unable to open $j ".$silentfiles[$j].".";

	my @counts_j = (0) x $nrmsbins;
	my @rms_j;
	my @energies_j;

	my ($scorecol,$rmscol) = (-1,-1);
	while (my $line=<SIL>) {
		chomp $line;
		next if ($line !~ /^SCORE:/);

		my @fields = split ' ', $line;
		if ($line =~ /description/) {
			foreach my $i (0..$#fields) {
				if ($fields[$i] eq $scoretag) { $scorecol=$i; }
				if ($fields[$i] eq "rms1" || $fields[$i] eq "rms") { $rmscol=$i; }
			}
		} else {
			if ($scorecol==-1) {
				die "unable to find score column";
			}
			if ($rmscol==-1) {
				die "unable to find rms column";
			}
			my $score_i = $fields[$scorecol];
			my $rms_i = $fields[$rmscol];

			next if ($score_i eq "nan" || $rms_i eq "nan" || $score_i eq "-nan" || $rms_i eq "-nan");
			if ($score_i>$SCORECAP) { $score_i = $SCORECAP; }

			push @energies_j, $score_i;
			push @rms_j, $rms_i;
			$sumvalidstructs ++;

			foreach my $i (0..$nrmsbins-1) {
				if ($rms_i >= $RMSBINS[$i] && $rms_i < $RMSBINS[$i+1]) {
					$counts_j[$i]++;
				}
			}
		}
	}

	push @all_counts, \@counts_j;
	push @all_energies, \@energies_j;
	push @all_rms, \@rms_j;

	close (SIL);
}

my @descrims;
foreach my $i (0..$nsilent-1) {
	# count bins
	my @usebins = (1) x $nrmsbins;
	foreach my $j (0..$nrmsbins-2) {
		my ($countlow,$counthi) = (0,0);
		foreach my $lowbins (0..$j) { $countlow += $all_counts[$i]->[$lowbins]; }
		foreach my $hibins ($j+1..$nrmsbins-1) { $counthi += $all_counts[$i]->[$hibins]; }
		if ($countlow < $MIN_BIN_SIZE || $counthi < $MIN_BIN_SIZE ) {
			$usebins[$j] = 0;
		}
	}

	my $nvalidbins = 0;
	foreach my $i (0..$nrmsbins-2) { $nvalidbins += $usebins[$i]; }


	# normalization
	my @sortenergies = sort {$a <=> $b} @{ $all_energies[$i] };
	my $scalefactor = $sortenergies[floor(0.5+0.95*scalar(@sortenergies))] - $sortenergies[floor(0.5+0.05*scalar(@sortenergies))];
	my $lowE = $sortenergies[0];

	if ($scalefactor == 0) { $scalefactor = 1; } # all scores > 0 cause this

	my $binavg=0;
	foreach my $j (0..$nrmsbins-2) {
		next if (!$usebins[$j]);

		my $boltz_num = 0;
		my $boltz_denom = 0;

		foreach my $k (0..scalar(@{$all_energies[$i]})-1) {
			my $boltzE = exp( -($all_energies[$i]->[$k] - $lowE) / ($TEMP * $scalefactor) );
			$boltz_denom += $boltzE;
			if ($all_rms[$i]->[$k] < $RMSBINS[$j+1]) {
				$boltz_num += $boltzE;
			}
		}

		my $binscore = ($boltz_num)/$boltz_denom;
		$binavg += $binscore;
	}
	if ($nvalidbins > 0) {
		my $descrim_i = ($binavg/$nvalidbins);
		$descrim_i = min($descrim_i , $cap);
		$descrim_i = max($descrim_i , -$cap);
		$avg_descrim += $descrim_i;
		push @descrims, $descrim_i;
	} else {
		push @descrims, 0;
	}
	$sumvalidbins += $nvalidbins;
}

print $avg_descrim/$nsilent." ".$sumvalidbins." ".$sumvalidstructs;
print "\n";
print "\n";

if ($verbose) {
	foreach my $i (0..$nsilent-1) {
		print $silentfiles[$i].' '.$descrims[$i]."\n";
	}
}

exit 0;

########
# subs #
########
sub min ($$) { $_[$_[0] > $_[1]] }
sub max ($$) { $_[$_[0] < $_[1]] }

