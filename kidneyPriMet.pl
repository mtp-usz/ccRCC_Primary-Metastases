#!/usr/bin/perl -w
# Author: Abdullah Kahraman
# Date: 20.01.2020

###############################################################################
###############################################################################
### Analyses Kidney Prim-Met Matched FMI Data.                              ###
###############################################################################
###############################################################################

use strict;
use warnings;
use Getopt::Long;
no autovivification;

my (
    # variable for parameters which are read in from commandline
    $help,
    $in,
    $mutStats,
   );

##############################################################################
### read all needed parameters from commandline ##############################

&GetOptions(
    "help!"      => \$help,     # print this help
    "in=s"       => \$in,       # e.g. data.tsv
    "mut!"       => \$mutStats, # print out mutations statistics
) or die "\nTry \"$0 -h\" for a complete list of options\n\n";

##############################################################################

# help
if ($help) {printHelp(); exit}

##############################################################################
### SETTINGS #################################################################
##############################################################################

##############################################################################
### SUBROUTINES ##############################################################
############################################################################## 

###############################################################################
sub printHelp {
###############################################################################
    # prints a help about the using and parameters of this scripts 
    # (execute if user types commandline parameter -h)
    # param:  no paramaters
    # return: no return value

    my (
	$usage,
	$sourceCode,
	@rows,
	$row,
	$option,
	$scriptInfo,
	$example,
       );

    $usage = "$0\n";


    print "\nUsage: " .  $usage . "\n";

    print "Valid options are:\n\n";
    open(MYSELF, "$0") or
      die "Cannot read source code file $0: $!\n";
    $sourceCode .= join "", <MYSELF>;
    close MYSELF;
    $sourceCode =~ s/^.+?\&GetOptions\(\n//s;
    $sourceCode =~ s/\n\).+$//s;
    @rows = split /\n/, $sourceCode;
    foreach $row (@rows){
        $option = $row;
	$option =~ s/\s+\"//g;
	$option =~ s/\"\s.+\#/\t\#/g;
	$option =~ s/=./\t<value> [required]/;
	$option =~ s/:./\t<value> [optional]/;
	$option =~ s/!/\t<non value> [optional]/;

	$row =~ s/^.*//;
	print "\t";
	printf("%-1s%-30s%-30s\n", "-",$option,$row);

    } # end of foreach $row (@rows)
    print "\n";
    print "Options may be abreviated, e.g. -h for --help\n\n";

    $example  = "$0";
}
##############################################################################

sub assessKeys {
    my ($hash1, $hash2) = @_;
    my (%hash3, %hash4);
    foreach my $key (sort keys %$hash1) {
	$key =~ s/ .*//;
	$hash3{$key} = "";
    }
    foreach my $key (sort keys %$hash2) {
	$key =~ s/ .*//;
	$hash4{$key} = "";
    }
    
    my $common = "";
    my $diff1  = "";
    my $diff2  = "";
    foreach my $key (sort keys %hash3) {
	if(exists $hash4{$key}) {
	    $common .= ",$key";
	} else {
	    $diff1 .= ",$key";
	}
    }
    foreach my $key (sort keys %hash4) {
	if(!exists $hash3{$key}) {
	    $diff2 .= ",$key";
	}
    }
    $diff1 =~ s/^,//;
    $diff2 =~ s/^,//;
    $common =~ s/^,//;
    return "$common\t$diff1\t$diff2";
}
##############################################################################

sub assessMuts {
    my ($hash1, $hash2) = @_;
    my (%hash3, %hash4);
    
    my $svCommon = "";
    my $svDiff1  = "";
    my $svDiff2  = "";
    my $snvCommon = "";
    my $snvDiff1  = "";
    my $snvDiff2  = "";
    my $indelCommon = "";
    my $indelDiff1  = "";
    my $indelDiff2  = "";
    my $stopCommon = "";
    my $stopDiff1  = "";
    my $stopDiff2  = "";
    my $spliceCommon = "";
    my $spliceDiff1  = "";
    my $spliceDiff2  = "";
    my $nonsenseCommon = "";
    my $nonsenseDiff1  = "";
    my $nonsenseDiff2  = "";

    foreach my $key (sort keys %$hash1) {
	$key =~ s/ \(.*//;
	$hash3{$key} = "";
    }
    foreach my $key (sort keys %$hash2) {
	$key =~ s/ \(.*//;
	$hash4{$key} = "";
    }

    foreach my $key (sort keys %hash3) {
	if(exists $hash4{$key}) {
	    $svCommon .= ",$key";
	    if($key =~ /fs/) {
		$indelCommon .= ",$key";
		$nonsenseCommon .= ",$key";
	    } elsif($key =~ /\*/) {
		$stopCommon .= ",$key";
		$nonsenseCommon .= ",$key";
	    } elsif($key =~ /splice/) {
		$spliceCommon .= ",$key";
		$nonsenseCommon .= ",$key";
	    } else {
		$snvCommon .= ",$key";
	    }
	} else {
	    $svDiff1 .= ",$key";
	    if($key =~ /fs/) {
		$indelDiff1 .= ",$key";
		$nonsenseDiff1 .= ",$key";
	    } elsif($key =~ /\*/) {
		$stopDiff1 .= ",$key";
		$nonsenseDiff1 .= ",$key";
	    } elsif($key =~ /splice/) {
		$spliceDiff1 .= ",$key";
		$nonsenseDiff1 .= ",$key";
	    } else {
		$snvDiff1 .= ",$key";
	    }
	}
    }
    foreach my $key (sort keys %hash4) {
	if(!exists $hash3{$key}) {
	    $svDiff2 .= ",$key";
	    if($key =~ /fs/) {
		$indelDiff2 .= ",$key";
		$nonsenseDiff2 .= ",$key";
	    } elsif($key =~ /\*/) {
		$stopDiff2 .= ",$key";
		$nonsenseDiff2 .= ",$key";
	    } elsif($key =~ /splice/) {
		$spliceDiff2 .= ",$key";
		$nonsenseDiff2 .= ",$key";
	    } else {
		$snvDiff2 .= ",$key";
	    }
	}
    }
    $svCommon =~ s/^,//;
    $svDiff1 =~ s/^,//;
    $svDiff2 =~ s/^,//;
    $snvCommon =~ s/^,//;
    $snvDiff1 =~ s/^,//;
    $snvDiff2 =~ s/^,//;
    $indelCommon =~ s/^,//;
    $indelDiff1 =~ s/^,//;
    $indelDiff2 =~ s/^,//;
    $stopCommon =~ s/^,//;
    $stopDiff1 =~ s/^,//;
    $stopDiff2 =~ s/^,//;
    $spliceCommon =~ s/^,//;
    $spliceDiff1 =~ s/^,//;
    $spliceDiff2 =~ s/^,//;
    $nonsenseCommon =~ s/^,//;
    $nonsenseDiff1 =~ s/^,//;
    $nonsenseDiff2 =~ s/^,//;
    return "$svCommon\t$svDiff1\t$svDiff2\t".
	"$snvCommon\t$snvDiff1\t$snvDiff2\t".
	"$indelCommon\t$indelDiff1\t$indelDiff2\t".
	"$stopCommon\t$stopDiff1\t$stopDiff2\t".
	"$spliceCommon\t$spliceDiff1\t$spliceDiff2\t".
	"$nonsenseCommon\t$nonsenseDiff1\t$nonsenseDiff2";
}
##############################################################################
### END OF SUBROUTINES########################################################
############################################################################## 


############
### MAIN ###
############

if(!defined $in) {
    print STDERR "\nPlease provide an input file. More information with $0 -help.\n\n";
    exit;
} else {
    open(F, $in) or die "ERROR: Failed to open input file: $in\n";
    my %data = ();
    my $header = <F>;
    while(my $l = <F>) {
	chomp($l);
	next if($l =~ /^#/);
	my ($uszId, $groupId, $externalId, $fmiId, $ngsNo, $dob, $gender, $disease, $specimen, $status,
	    $sv, $cna, $fusion, $svVUS, $cnaVUS, $fusionVUS,
	    $tmb, $msStatus, $msScore, $lohStatus, $loh, $dna, $tcUSZ, $tcFMI, $tcFMIstatus) = split(/\t/, $l);
	# skip following cases
	next if(
	    $uszId =~ /LC2019\.224\./ or
	    $uszId =~ /LC2019\.259\./ or
	    $uszId =~ /LC2019\.261\./ or
	    $uszId =~ /LC2019\.268\./ or
	    $uszId =~ /LC2019\.269\./ or
	    $uszId =~ /LC2019\.282\./ or
	    $uszId =~ /LC2019\.288\./ or		
	    $uszId =~ /LC2019\.289\./ or
	    $uszId =~ /LC2020\.504\./);
	if($uszId =~ /LC20\d\d\.(\d+)\.1\.(\d+).1/) {
	    my $uszNumber = $1;
	    my $uszBlock  = $2;
	    if($fmiId ne "") {
		if($specimen eq "Kidney") {
		    if(!exists $data{$uszNumber}->{"Primary"}) {
			$data{$uszNumber}->{"Primary"} = "$specimen\t$dob\t$gender\t$disease\t$sv\t$cna\t$fusion\t$svVUS\t$cnaVUS\t$fusionVUS\t$tmb\t$msScore\t$loh";
		    } else {
			print STDERR "WARNING: Multiple primary: $l\n";
		    }
		} else {
		    # use last metastases
		    $data{$uszNumber}->{"Metastase"} = "$specimen\t$dob\t$gender\t$disease\t$sv\t$cna\t$fusion\t$svVUS\t$cnaVUS\t$fusionVUS\t$tmb\t$msScore\t$loh";
		}
	    }
	} else {
	    die "ERROR: $uszId has wrong format\n";
	}	
    }

    if(!defined $mutStats) {
	print "#USZid\tSpecimen\tDOB\tGender\tDisease\t".
	    "SVcommon\tSVdiffP\tSVdiffM\t".
	    "CNAcommon\tCNAdiffP\tCNAdiffM\t".
	    "fusionCommon\tfusionDiffP\tfusionDiffM\t".
	    "vusSVcommon\tvusSVdiffP\tvusSVdiffM\t".
	    "vusCNAcommon\tvusCNAdiffP\tvusCNAdiffM\t".
	    "VUSfusionCommon\tVUSfusionDiffP\tVUSfusionDiffM\t".
	    "AllCommon\tAllDiffP\tAllDiffM\t".
	    "MSIp\tMSIm\tTMBp\tTMBm\tLOHp\tLOHm\n";
    } else {
	print "#USZid\tSpecimen\tDOB\tGender\tDisease\t".
	    "SVcommon\tSVdiffP\tSVdiffM\t".
	    "SNVcommon\tSNVdiffP\tSNVdiffM\t".
	    "indelCommon\tindelDiffP\tindelDiffM\t".
	    "stopCommon\tstopDiffP\tstopDiffM\t".
	    "spliceCommon\tspliceDiffP\tspliceDiffM\t".
	    "nonsenseCommon\tnonsenseDiffP\tnonsenseDiffM\t".
	    "vusSVcommon\tvusSVdiffP\tvusSVdiffM\t".
	    "vusSNVcommon\tvusSNVdiffP\tvusSNVdiffM\t".
	    "VUSindelCommon\tVUSindelDiffP\tVUSindelDiffM\t".
	    "VUSstopCommon\tVUSstopDiffP\tstopDiffM\t".
	    "VUSspliceCommon\tVUSspliceDiffP\tVUSspliceDiffM\t".
	    "VUSnonsenseCommon\tVUSnonsenseDiffP\tVUSnonsenseDiffM\n";

    }
    foreach my $uszNumber (sort keys %data) {
	next if(!exists $data{$uszNumber}->{"Metastase"} or !exists $data{$uszNumber}->{"Primary"});
	my $primary = $data{$uszNumber}->{"Primary"};
	my $metastase = $data{$uszNumber}->{"Metastase"};
	my $specimenP = "";
	my $dobP = "";
	my $genderP = "";
	my $diseaseP = "";
	my $svP = "";
	my $cnaP = "";
	my $fusionP = "";
	my $svVUSp = "";
	my $cnaVUSp = "";
	my $fusionVUSp = "";
	my $tmbP = "";
	my $msScoreP = "";
	my $lohP = "";
	my $specimenM = "";
	my $dobM = "";
	my $genderM = "";
	my $diseaseM = "";
	my $svM = "";
	my $cnaM = "";
	my $fusionM = "";
	my $svVUSm = "";
	my $cnaVUSm = "";
	my $fusionVUSm = "";
	my $tmbM = "";
	my $msScoreM = "";
	my $lohM = "";
	($specimenP, $dobP, $genderP, $diseaseP, $svP, $cnaP, $fusionP, $svVUSp, $cnaVUSp, $fusionVUSp, $tmbP, $msScoreP, $lohP) = split(/\t/, $primary)
	    if(defined $primary); 
	($specimenM, $dobM, $genderM, $diseaseM, $svM, $cnaM, $fusionM, $svVUSm, $cnaVUSm, $fusionVUSm, $tmbM, $msScoreM, $lohM) = split(/\t/, $metastase)
	    if(defined $metastase);
	    
	my %svP = map{ $_ => ""} split(/,/, $svP);
	my %svM = map{ $_ => ""} split(/,/, $svM);
	my %cnaP = map{ $_ => ""} split(/,/, $cnaP);
	my %cnaM = map{ $_ => ""} split(/,/, $cnaM);
	my %fusionP = map{ $_ => ""} split(/,/, $fusionP);
	my %fusionM = map{ $_ => ""} split(/,/, $fusionM);
	my %svVUSp = map{ $_ => ""} split(/,/, $svVUSp);
	my %svVUSm = map{ $_ => ""} split(/,/, $svVUSm);
	my %cnaVUSp = map{ $_ => ""} split(/,/, $cnaVUSp);
	my %cnaVUSm = map{ $_ => ""} split(/,/, $cnaVUSm);
	my %fusionVUSp = map{ $_ => ""} split(/,/, $fusionVUSp);
	my %fusionVUSm = map{ $_ => ""} split(/,/, $fusionVUSm);
	my $p = "$svP,$cnaP,$fusionP,$svVUSp,$cnaVUSp,$fusionVUSp";
	my $m = "$svM,$cnaM,$fusionM,$svVUSm,$cnaVUSm,$fusionVUSm";
	$p =~ s/^,+//;
	$p =~ s/,+$//;
	$p =~ s/,+/,/g;
	$m =~ s/^,+//;
	$m =~ s/,+$//;
	$m =~ s/,+/,/g;
	my %p = map{ $_ => ""} split(/,/, $p);
	my %m = map{ $_ => ""} split(/,/, $m);
	if(($dobP ne $dobM and $dobP ne "" and $dobM ne "") or
	   ($genderP ne $genderM and $genderP ne "" and $genderM ne "") or
	   ($diseaseP ne $diseaseM and $diseaseP eq "" and $diseaseM eq "")) {
	    print "Discordance\n";
	} else {
	    $specimenM = $specimenP if($specimenM eq "");
	    $dobM = $dobP if($dobM eq "");
	    $genderM = $genderP if($genderM eq "");
	    $diseaseM = $diseaseP if($diseaseM eq "");
	    print "LC2019.$uszNumber\t$specimenM\t$dobP\t$genderP\t$diseaseP\t";
	    if(!defined $mutStats) {
		print &assessKeys(\%svP, \%svM)."\t".&assessKeys(\%cnaP, \%cnaM)."\t".&assessKeys(\%fusionP, \%fusionM)."\t".
		&assessKeys(\%svVUSp, \%svVUSm)."\t".&assessKeys(\%cnaVUSp, \%cnaVUSm)."\t".&assessKeys(\%fusionVUSp, \%fusionVUSm)."\t".
		&assessKeys(\%p, \%m).
		"\t$msScoreP\t$msScoreM\t$tmbP\t$tmbM\t$lohP\t$lohM\n";
	    } else {
		print &assessMuts(\%svP, \%svM)."\t".&assessMuts(\%svVUSp, \%svVUSm)."\n";
	    }
	    print STDERR $primary if(!defined $lohP);
	}
    }
}
