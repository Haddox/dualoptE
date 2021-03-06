#!/usr/bin/perl
##
##
###############################################################################

use strict;
use POSIX;
use Carp;
use Getopt::Long;
use constant TINY => 1e-16;


my ($median) = 0;
GetOptions ("median=i" => \$median);

# constants
my ($ALPHA,$BETA,$GAMMA)=(1.0,0.5,2.0);
my $STRUCTCOUNTER=1;

# PARAMS
my $verbose = 0;
my $kT = 0.5;

if ($#ARGV<0) {
	print STDERR "usage: $0 <perres_fasols>\n";
	exit -1;
}

my $infile = shift @ARGV;

# E_transfer (vapor->water)
my %dgrefs = (
	'A' => 1.94,
	'C' => -1.24,
	'D' => -10.95,
	'E' => -10.24,
	'F' => -0.76,
#	'G' => 2.39,
	'H' => -10.27,
	'I' => 2.15,
	'K' => -9.52,
	'L' => 2.28,
	'M' => -1.48,
	'N' => -9.68,
	'Q' => -9.38,
	'R' => -19.92,
	'S' => -5.06,
	'T' => -4.88,
	'V' => 1.99,
	'W' => -5.88,
	'Y' => -6.1
);

# SA
my %SAs = (
	'L' => 121.2,
	'I' => 129.6,
	'V' => 104.8,
	'F' => 176,
	'M' => 144.8,
	'W' => 219.2,
	'A' => 41.2,
	'C' => 82.4,
#	'G' => 0,
	'Y' => 187.2,
	'T' => 89.6,
	'S' => 56.8,
	'H' => 135.6,
	'Q' => 125.2,
	'K' => 162.8,
	'N' => 93.6,
	'E' => 130.4,
	'D' => 90,
	'R' => 198
);

my (%dgros,%ljros);
my (%solcounts,%vdwcounts,%burcounts);

open (FASOL, $infile) || die "Unable to open $infile.";
my @fasols = <FASOL>;
chomp (@fasols);

foreach my $line ( @fasols ) {
	my @fields = split ' ', $line;
	my $aa = $fields[4];
	my $Esol = $fields[5];
	my $Elj = $fields[6];
	my $burial = $fields[7];

	if (!defined $solcounts{$aa}) { $solcounts{$aa} = []; }
	if (!defined $vdwcounts{$aa}) { $vdwcounts{$aa} = []; }
	if (!defined $burcounts{$aa}) { $burcounts{$aa} = []; }

	if ($burial eq "nan" || $burial eq "-nan") { next; }

	push @{$solcounts{$aa}}, $Esol;
	push @{$vdwcounts{$aa}}, $Elj;
	push @{$burcounts{$aa}}, $burial;
}

my ($xysum,$xxsum,$yysum,$xsum,$ysum) = (0,0,0);

# least fit line
foreach my $aa (keys %dgrefs) {
	my $nelts = scalar(@{$burcounts{$aa}});
	my @idxes =  sort { $burcounts{$aa}->[$b] <=> $burcounts{$aa}->[$a] } 0..$nelts-1;

	my @sortedEsol = @{$solcounts{$aa}}[ @idxes ];
	my @sortedLJ = @{$vdwcounts{$aa}}[ @idxes ];

	my ($sumLK,$sumLJ) = (0,0);
	if ($median) {
		my $start = int(($nelts-$median)/2);
		my $stop = $start+$median-1;
		if ($start < 0) { $start = 0; }
		if ($stop > $nelts-1) { $stop = $nelts-1; }
		foreach my $i ($start..$stop) {
			$sumLK += -$sortedEsol[$i];
			$sumLJ += $sortedLJ[$i];
		}
		$sumLK /= ($stop-$start+1);
		$sumLJ /= ($stop-$start+1);
	} else {
		my $N = 5; ## $nelts/5;
		if ($nelts < $N) { $N = $nelts; }
	
		foreach my $i (0..$N-1) {
			$sumLK += -$sortedEsol[$i]; # - $sortedLJ[$i];
			$sumLJ += $sortedLJ[$i];
		}
		$sumLK /= $N;
		$sumLJ /= $N;
	}
	
	$ljros{ $aa } = $sumLJ;
	$dgros{ $aa } = $sumLK;

	my $x = $dgrefs{ $aa } / $SAs{$aa};
	my $y = $dgros{ $aa } / $SAs{$aa};

	$xsum += $x;
	$ysum += $y;
	$xysum += $x*$y;
	$xxsum += $x*$x;
	$yysum += $y*$y;
}

$xsum /= scalar(keys %dgrefs);
$ysum /= scalar(keys %dgrefs);
$xysum /= scalar(keys %dgrefs);
$xxsum /= scalar(keys %dgrefs);
$yysum /= scalar(keys %dgrefs);

my $M = ($xysum - $xsum*$ysum) / ($yysum - $ysum*$ysum);
my $B = $xsum - $M*$ysum;

my ($linLoss,$absLoss,$denom) = (0,0,0);
foreach my $aa (keys %dgrefs) {
	# dgref/SA = (eps)*dgros/SA - C*(1-1/eps)
	my $dgrefsa = $dgrefs{$aa}/$SAs{$aa};
	my $dgrefsa_fit = $dgros{$aa}/$SAs{$aa};
	my $linErr = ( $dgrefsa  - ($M * $dgrefsa_fit + $B));
	my $absErr = ( $dgrefsa  - $dgrefsa_fit );
	$linLoss += $linErr * $linErr;
	$absLoss += $absErr * $absErr;
	$denom += $dgrefsa * $dgrefsa;
}
$linLoss /= $denom;
$absLoss /= $denom;

# back calculate epsilon and c
my $epsilon = $M;
my $C = -$B / (1-1/$M);

print $linLoss." ".$absLoss." ".scalar(@fasols)."\n";
print "eps/c = ".$epsilon." ".$C."\n";
print "M/B = ".$M." ".$B."\n";
foreach my $aa (keys %dgrefs) {
	print $aa." ".$dgros{ $aa }/$SAs{$aa}." ".$dgrefs{$aa}/$SAs{$aa}." "." ".scalar(@{$burcounts{$aa}})."\n";
}