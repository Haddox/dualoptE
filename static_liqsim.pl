#!/usr/bin/perl
##
##
###############################################################################

use strict;
use POSIX;

if ($#ARGV < 0) {
	print STDERR "usage: $0 <fracN> <scorefile>\n";
	exit -1;
}

my $FRACN = shift @ARGV;

my $TTEMP = 30.0;  #from hahnbeom
my $TVOL = 15.0;   #from hahnbeom
my $CTEMP = 1.0/$TTEMP;
my $CVOL = 1.0/$TVOL;
my $P0 = 1;
my $N = 64;

my %ref_temp=(
	'ethane' => 185,
	'propane' => 231,
	'isobutane' => 298,
	'benzene' => 298,
	'toluene' => 298,
	'methanol' => 298,
	'ethanol' => 298,
	'phenol' => 298,
	'CH3SH' => 279,
	'CH3SCH3' => 298,
	'CH3CH2SH' => 298,
	'acetamide' => 373,
	'DME' => 248,
	'NMF' => 298,
	'ethanal' => 298,
	'propanone' => 298,
	't-BuOH' => 298
);

my %ref_dens=(
	'ethane' => 0.546,
	'propane' => 0.581,
	'isobutane' => 0.551,
	'benzene' => 0.880,
	'toluene' => 0.864,
	'methanol' => 0.786,
	'ethanol' => 0.785,
	'phenol' => 1.058,
	'CH3SH' => 0.888,
	'CH3SCH3' => 0.842,
	'CH3CH2SH' => 0.833,
	'acetamide' => 0.981,
	'DME' => 0.735,
	'NMF' => 1.011,
	'ethanal' => 0.772,
	'propanone' => 0.784,
	't-BuOH' => 0.781
);

my %ref_hvap=(
	'ethane' => 3.520,
	'propane' => 4.490,
	'isobutane' => 4.570,
	'benzene' => 8.080,
	'toluene' => 9.080,
	'methanol' => 8.950,
	'ethanol' => 10.110,
	'phenol' => 13.820,
	'CH3SH' => 5.870,
	'CH3SCH3' => 6.610,
	'CH3CH2SH' => 6.580,
	'acetamide' => 14.200,
	'DME' => 5.140,
	'NMF' => 13.430,
	'ethanal' => 6.240,
	'propanone' => 7.480,
	't-BuOH' => 11.140
);

##
## computing predicted Hsim
##   Hvap_pred = (Hvap_static - Hvap_static_ref) * hvap_scale + Hvap_sim_ref
my %hvap_sim_ref = (
	'ethane' => 3.472,
	'propane' => 4.249,
	'isobutane' => 4.226,
	'benzene' => 7.709,
	'toluene' => 9.546,
	'methanol' => 8.240,
	'ethanol' => 9.044,
	'phenol' => 14.438,
	'CH3SH' => 4.770,
	'CH3SCH3' => 5.896,
	'CH3CH2SH' => 4.130,
	'acetamide' => 12.144,
	'DME' => 4.129,
	'NMF' => 10.807,
	'ethanal' => 4.261,
	'propanone' => 6.454,
	't-BuOH' => 11.281
);

## UPDATE!
my %hvap_static_ref = (
	'benzene' => 7.55020167569161,
	'propanone' => 7.27490097842273,
	'methanol' => 9.85617239918778,
	't-BuOH' => 12.0574296888866,
	'acetamide' => 15.2184695528567,
	'isobutane' => 4.93379154147492,
	'DME' => 4.37436952762104,
	'CH3CH2SH' => 5.70448244270257,
	'toluene' => 9.4087692878096,
	'CH3SH' => 4.96758337550354,
	'phenol' => 15.7498646224601,
	'NMF' => 12.0379584675997,
	'CH3SCH3' => 6.29706902302599,
	'ethane' => 3.90301944890122,
	'ethanol' => 10.9262533579892,
	'ethanal' => 4.87559570342887,
	'propane' => 4.80979422243894,
);

##
## computing predicted Dens
##   dens_pred = (Hvap_static - Hvap_static_ref) * hvap_scale + Hvap_sim_ref
my %dens_sim_ref = (
	'ethane' => 0.584,
	'propane' => 0.624,
	'isobutane' => 0.585,
	'benzene' => 0.921,
	'toluene' => 0.922,
	'methanol' => 0.717,
	'ethanol' => 0.756,
	'phenol' => 1.126,
	'CH3SH' => 0.960,
	'CH3SCH3' => 0.870,
	'CH3CH2SH' => 0.814,
	'acetamide' => 1.003,
	'DME' => 0.728,
	'NMF' => 0.981,
	'ethanal' => 0.744,
	'propanone' => 0.838,
	't-BuOH' => 0.810
);

my %dens_static_ref = (
	'benzene' => 0.93236507936508,
	'propanone' => 0.897307317073171,
	'methanol' => 0.803556097560976,
	't-BuOH' => 0.869597883597884,
	'acetamide' => 0.962029268292683,
	'isobutane' => 0.661275862068966,
	'DME' => 0.796149425287356,
	'CH3CH2SH' => 0.91224358974359,
	'toluene' => 0.937285714285714,
	'CH3SH' => 0.950949367088608,
	'phenol' => 1.11842857142857,
	'NMF' => 1.03906329113924,
	'CH3SCH3' => 0.914992063492063,
	'ethane' => 0.648711711711712,
	'ethanol' => 0.838303191489362,
	'ethanal' => 0.812468292682928,
	'propane' => 0.686676056338028
);

##
## relative weights
my %weights = (
	'ethane' => 1,
	'propane' => 1,
	'isobutane' => 1,
	'benzene' => 1,
	'toluene' => 1,
	'methanol' => 1,
	'ethanol' => 1,
	'phenol' => 1,
	'CH3SH' => 1,
	'CH3SCH3' => 1,
	'CH3CH2SH' => 1,
	'acetamide' => 1,
	'DME' => 0.5,
	'NMF' => 1,
	'ethanal' => 0.5,
	'propanone' => 0.5,
	't-BuOH' => 1
);

my $dens_wt = 0.5;

## ???
my $hvap_scale = 2.5;
my $dens_scale = 4.0;

my %calc_dens;
my %calc_hvap;

my $nstruct=0;

my (%allDs, %allVs, %allEs);
my (%Eref, %Vref);

foreach my $scfile (@ARGV) {
	open (SIL, $scfile) || die "Unable to open ".$scfile.".";

	my @scsplit = split(/\./, $scfile);
	my $tag = $scsplit[0];
	$tag =~ s/.*\///g;

	if (!defined $allDs{$tag}) {
		$allDs{$tag} = [];
		$allVs{$tag} = [];
		$Eref{$tag} = 1e6;
		$Vref{$tag} = 1e6;
	}

	my ($scorecol,$volcol,$denscol,$corrcol,$intraAcol,$intraRcol) = (-1,-1,-1,-1,-1);

	while (my $line=<SIL>) {
		chomp $line;
		next if ($line !~ /^SCORE:/);

		my @fields = split ' ', $line;
		if ($line =~ /description/) {
			foreach my $i (0..$#fields) {
				if ($fields[$i] eq "score") { $scorecol=$i; }
				if ($fields[$i] eq "fa_intra_atr_xover4") { $intraAcol=$i; }
				if ($fields[$i] eq "fa_intra_rep_xover4") { $intraRcol=$i; }
				if ($fields[$i] eq "LJcorr") { $corrcol=$i; }
				if ($fields[$i] eq "volume") { $volcol=$i; }
				if ($fields[$i] eq "density") { $denscol=$i; }
			}
		} else {
			if ($scorecol==-1) {
				die "unable to find score column";
			}
			if ($volcol==-1) {
				die "unable to find volume column";
			}
			if ($denscol==-1) {
				die "unable to find dens column";
			}

			my $score_i = ($fields[$scorecol] + $fields[$corrcol] - $fields[$intraAcol] - $fields[$intraRcol])*0.5;
			my $vol_i = $fields[$volcol];
			my $dens_i = $fields[$denscol];

			#if ($Vref>$vol_i) { $Vref=$vol_i; }
			if ($Eref{$tag}>$score_i) {
				$Eref{$tag}=$score_i;
				$Vref{$tag}=$vol_i;
			}

			push @{$allDs{$tag}}, $dens_i;
			push @{$allVs{$tag}}, $vol_i;
			push @{$allEs{$tag}}, $score_i;
		}
	}
	close (SIL);
}

foreach my $tag (keys %allDs) {
	my @allD = @{$allDs{$tag}};
	my @allV = @{$allVs{$tag}};
	my @allE = @{$allEs{$tag}};

	my $Eref = $Eref{$tag};
	my $Vref = $Vref{$tag};

	my $TEMP = $ref_temp{$tag};

	my ($Psum)=(0);
	my (@allP, @allLogP);

	my $TOPN = ceil($FRACN*scalar(@allD));

	foreach my $i (0..$#allD) {
		my $score_i = $allE[$i];
		my $vol_i = $allV[$i];
		my $dens_i = $allD[$i];

		my $Eaft = ( $score_i );

		my $del_volume = ($vol_i - $Vref);

		my $beta = 1.0/(($TEMP+0.00100)*0.0019872065);
		my $convert = 0.000014586;
		my $arg1 = $beta*($Eaft-$Eref) + $P0*$del_volume*$convert;
		my $arg2 = -($N) * log( $vol_i/$Vref );

		my $logP = $CTEMP*$arg1 + $CVOL*$arg2;
		my $P = exp( -($logP) );

		#print " ".$logP." ".$dens_i."\n";
		$Psum += $P;
		push @allP, $P;
		push @allLogP, $logP;
	}

	my @sortLogPidx = sort { $allLogP[$a] <=> $allLogP[$b] } 0..$#allLogP;
	my @sortLogP = @allLogP[ @sortLogPidx ];
	#print STDERR " ".$sortLogP[0]." , ".$sortLogP[1]." , ".$sortLogP[2]." , ... , ".$sortLogP[$#allLogP]."\n";
	my @sortE = @allE[ @sortLogPidx ];
	my @sortD = @allD[ @sortLogPidx ];
	my @sortP = @allP[ @sortLogPidx ];

	#print STDERR "Have ".scalar(@allP)."\n";
	my ($Eavg,$Davg)=(0,0,0);
	foreach my $i (0..($TOPN-1)) {
		#print " ".$sortLogP[$i]." ".($sortP[$i]/$Psum)." ".$sortD[$i]." ".$sortE[$i]."\n";
		$Eavg += $sortE[$i]/$TOPN;
		$Davg += $sortD[$i]/$TOPN;
	}

	my @sortTopE = sort @sortE[0 .. ($TOPN-1)];
	my @sortTopD = sort @sortD[0 .. ($TOPN-1)];
	my $Emed = $sortTopE[$TOPN/2];
	my $Dmed = $sortTopD[$TOPN/2];


	my ($Ewt,$Vwt,$Dwt)=(0,0,0);
	foreach my $i (0..$#allP) {
		$Ewt += $allE[$i]*$allP[$i]/$Psum;
		$Dwt += $allD[$i]*$allP[$i]/$Psum;
	}

	#
	my $Erpt = -$Ewt/$N + 0.0019872965*$TEMP;
	my $Drpt = $Davg;

	$calc_hvap{ $tag } = $Erpt;
	$calc_dens{ $tag } = $Drpt;

	$nstruct += scalar(@allD);
}

# compute error
# a)
my $errH = 0;
my $errD = 0;
my $wtSUM = 0;
my (%hvap_pred,%dens_pred);
foreach my $tag (keys %calc_hvap) {

	##
	## NOTE: BASELINE FOR HVAP and DENS ARE DIFFERENT!
	$hvap_pred{ $tag } = ($calc_hvap{ $tag } - $hvap_static_ref{$tag}) * $hvap_scale  + 0.5*( $hvap_static_ref{$tag} + $hvap_sim_ref{$tag} );
	$dens_pred{ $tag } = ($calc_dens{ $tag } - $dens_static_ref{$tag}) * $dens_scale  + 0.5*( $dens_static_ref{$tag} + $dens_sim_ref{$tag} );

	my $delH = ($hvap_pred{ $tag }-$ref_hvap{ $tag });
	my $delD = ($dens_pred{ $tag }-$ref_dens{ $tag });

	$errH += $weights{$tag} * abs($delH) / $ref_hvap{ $tag };
	$errD += $weights{$tag} * abs($delD) / $ref_dens{ $tag };

	$wtSUM += $weights{$tag};
}
$errH = ( $errH / $wtSUM );
$errD = ( $errD / $wtSUM );

my $error = $dens_wt*$errD+$errH;

print "$error $nstruct\n";
foreach my $tag (keys %calc_hvap) {
	print "$tag   ".
		$hvap_pred{ $tag }." ".$calc_hvap{ $tag }." ".$ref_hvap{ $tag }."   ".
		$dens_pred{ $tag }." ".$calc_dens{ $tag }." ".$ref_dens{ $tag }."\n";
}


exit 0;

########
# subs #
########
sub min ($$) { $_[$_[0] > $_[1]] }
sub max ($$) { $_[$_[0] < $_[1]] }

