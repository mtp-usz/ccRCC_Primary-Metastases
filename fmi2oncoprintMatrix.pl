#!/usr/bin/perl -w
# Author: Abdullah Kahraman
# Date: 17.04.2021

###############################################################################
###############################################################################
### Takes FMI data and converts it to an OncoPrint Matrix Input data.       ###
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
    $vus,
    $snvSplit,
    $minPerc,
    );

##############################################################################
### read all needed parameters from commandline ##############################

&GetOptions(
    "help!"      => \$help,     # print this help
    "in=s"       => \$in,       # input file, e.g. kidney_priMet_dual.tsv
    "vus!"       => \$vus,      # include also VUS
    "snv!"       => \$snvSplit, # split SNVs into Missense, Nonsense, Frameshift, Splice site and Promoter
    "min=i"      => \$minPerc,  # minimum percentage of mutations per sample, e.g. 10
) or die "\nTry \"$0 -h\" for a complete list of options\n\n";

##############################################################################

# help
if ($help) {printHelp(); exit}

##############################################################################
### SETTINGS #################################################################
##############################################################################

my %fmiGenes =
   ("ABL1" => "ABL1",
    "ACVR1B" => "ACVR1B",
    "AKT1" => "AKT1",
    "AKT2" => "AKT2",
    "AKT3" => "AKT3",
    "ALK" => "ALK",
    "ALOX12B" => "ALOX12B",
    "APC" => "APC",
    "AR" => "AR",
    "ARAF" => "ARAF",
    "ARFRP1" => "ARFRP1",
    "ARID1A" => "ARID1A",
    "ASXL1" => "ASXL1",
    "ATM" => "ATM",
    "ATR" => "ATR",
    "ATRX" => "ATRX",
    "AURKA" => "AURKA",
    "AURKB" => "AURKB",
    "AXIN1" => "AXIN1",
    "AXL" => "AXL",
    "BAP1" => "BAP1",
    "BARD1" => "BARD1",
    "BCL2" => "BCL2",
    "BCL2L1" => "BCL2L1",
    "BCL2L2" => "BCL2L2",
    "BCL6" => "BCL6",
    "BCOR" => "BCOR",
    "BCORL1" => "BCORL1",
    "BCR" => "BCR",
    "BRAF" => "BRAF",
    "BRCA1" => "BRCA1",
    "BRCA2" => "BRCA2",
    "BRD4" => "BRD4",
    "BRIP1" => "BRIP1",
    "BTG1" => "BTG1",
    "BTG2" => "BTG2",
    "BTK" => "BTK",
    "C11orf30" => "EMSY",
    "EMSY" => "EMSY",
    "C17orf39" => "GID4",
    "GID4" => "GID4",
    "CALR" => "CALR",
    "CARD11" => "CARD11",
    "CASP8" => "CASP8",
    "CBFB" => "CBFB",
    "CBL" => "CBL",
    "CCND1" => "CCND1",
    "CCND2" => "CCND2",
    "CCND3" => "CCND3",
    "CCNE1" => "CCNE1",
    "CD22" => "CD22",
    "CD70" => "CD70",
    "CD74" => "CD74",
    "CD79A" => "CD79A",
    "CD79B" => "CD79B",
    "CDC73" => "CDC73",
    "CDH1" => "CDH1",
    "CDK12" => "CDK12",
    "CDK4" => "CDK4",
    "CDK6" => "CDK6",
    "CDK8" => "CDK8",
    "CDKN1A" => "CDKN1A",
    "CDKN1B" => "CDKN1B",
    "CDKN2A" => "CDKN2A",
    "CDKN2B" => "CDKN2B",
    "CDKN2C" => "CDKN2C",
    "CEBPA" => "CEBPA",
    "CHEK1" => "CHEK1",
    "CHEK2" => "CHEK2",
    "CIC" => "CIC",
    "CREBBP" => "CREBBP",
    "CRKL" => "CRKL",
    "CSF1R" => "CSF1R",
    "CSF3R" => "CSF3R",
    "CTCF" => "CTCF",
    "CTNNA1" => "CTNNA1",
    "CTNNB1" => "CTNNB1",
    "CUL3" => "CUL3",
    "CUL4A" => "CUL4A",
    "CXCR4" => "CXCR4",
    "CYP17A1" => "CYP17A1",
    "DAXX" => "DAXX",
    "DDR1" => "DDR1",
    "DDR2" => "DDR2",
    "DIS3" => "DIS3",
    "DNMT3A" => "DNMT3A",
    "DOT1L" => "DOT1L",
    "EED" => "EED",
    "EGFR" => "EGFR",
    "EP300" => "EP300",
    "EPHA3" => "EPHA3",
    "EPHB1" => "EPHB1",
    "EPHB4" => "EPHB4",
    "ERBB2" => "ERBB2",
    "ERBB3" => "ERBB3",
    "ERBB4" => "ERBB4",
    "ERCC4" => "ERCC4",
    "ERG" => "ERG",
    "ERRFI1" => "ERRFI1",
    "ESR1" => "ESR1",
    "ETV4" => "ETV4",
    "ETV5" => "ETV5",
    "ETV6" => "ETV6",
    "EWSR1" => "EWSR1",
    "EZH2" => "EZH2",
    "EZR" => "EZR",
    "FAM123B" => "AMER1",
    "AMER1" => "AMER1",
    "FAM46C" => "FAM46C",
    "FANCA" => "FANCA",
    "FANCC" => "FANCC",
    "FANCG" => "FANCG",
    "FANCL" => "FANCL",
    "FAS" => "FAS",
    "FBXW7" => "FBXW7",
    "FGF10" => "FGF10",
    "FGF12" => "FGF12",
    "FGF14" => "FGF14",
    "FGF19" => "FGF19",
    "FGF23" => "FGF23",
    "FGF3" => "FGF3",
    "FGF4" => "FGF4",
    "FGF6" => "FGF6",
    "FGFR1" => "FGFR1",
    "FGFR2" => "FGFR2",
    "FGFR3" => "FGFR3",
    "FGFR4" => "FGFR4",
    "FH" => "FH",
    "FLCN" => "FLCN",
    "FLT1" => "FLT1",
    "FLT3" => "FLT3",
    "FOXL2" => "FOXL2",
    "FUBP1" => "FUBP1",
    "GABRA6" => "GABRA6",
    "GATA3" => "GATA3",
    "GATA4" => "GATA4",
    "GATA6" => "GATA6",
    "GNA11" => "GNA11",
    "GNA13" => "GNA13",
    "GNAQ" => "GNAQ",
    "GNAS" => "GNAS",
    "GRM3" => "GRM3",
    "GSK3B" => "GSK3B",
    "H3F3A" => "H3F3A",
    "HDAC1" => "HDAC1",
    "HGF" => "HGF",
    "HNF1A" => "HNF1A",
    "HRAS" => "HRAS",
    "HSD3B1" => "HSD3B1",
    "ID3" => "ID3",
    "IDH1" => "IDH1",
    "IDH2" => "IDH2",
    "IGF1R" => "IGF1R",
    "IKBKE" => "IKBKE",
    "IKZF1" => "IKZF1",
    "INPP4B" => "INPP4B",
    "IRF2" => "IRF2",
    "IRF4" => "IRF4",
    "IRS2" => "IRS2",
    "JAK1" => "JAK1",
    "JAK2" => "JAK2",
    "JAK3" => "JAK3",
    "JUN" => "JUN",
    "KDM5A" => "KDM5A",
    "KDM5C" => "KDM5C",
    "KDM6A" => "KDM6A",
    "KDR" => "KDR",
    "KEAP1" => "KEAP1",
    "KEL" => "KEL",
    "KIT" => "KIT",
    "KLHL6" => "KLHL6",
    "KRAS" => "KRAS",
    "LTK" => "LTK",
    "LYN" => "LYN",
    "MAF" => "MAF",
    "MAP2K4" => "MAP2K4",
    "MAP3K1" => "MAP3K1",
    "MAP3K13" => "MAP3K13",
    "MAPK1" => "MAPK1",
    "MCL1" => "MCL1",
    "MDM2" => "MDM2",
    "MDM4" => "MDM4",
    "MED12" => "MED12",
    "MEF2B" => "MEF2B",
    "MEK1" => "MAP2K1",
    "MAP2K1" => "MAP2K1",
    "MEK2" => "MAP2K2",
    "MAP2K2" => "MAP2K2",
    "MEN1" => "MEN1",
    "MERTK" => "MERTK",
    "MET" => "MET",
    "MITF" => "MITF",
    "MKNK1" => "MKNK1",
    "MLH1" => "MLH1",
    "MLL" => "KMT2A",
    "KMT2A" => "KMT2A",
    "MLL2" => "KMT2D",
    "KMT2D" => "KMT2D",
    "MPL" => "MPL",
    "MRE11A" => "MRE11A",
    "MSH2" => "MSH2",
    "MSH3" => "MSH3",
    "MSH6" => "MSH6",
    "MST1R" => "MST1R",
    "MTAP" => "MTAP",
    "MTOR" => "MTOR",
    "MUTYH" => "MUTYH",
    "MYB" => "MYB",
    "MYC" => "MYC",
    "MYCL1" => "MYCL",
    "MYCL" => "MYCL",
    "MYCN" => "MYCN",
    "MYD88" => "MYD88",
    "NBN" => "NBN",
    "NF1" => "NF1",
    "NF2" => "NF2",
    "NFE2L2" => "NFE2L2",
    "NFKBIA" => "NFKBIA",
    "NKX2-1" => "NKX2-1",
    "NOTCH1" => "NOTCH1",
    "NOTCH2" => "NOTCH2",
    "NOTCH3" => "NOTCH3",
    "NPM1" => "NPM1",
    "NRAS" => "NRAS",
    "NT5C2" => "NT5C2",
    "NTRK1" => "NTRK1",
    "NTRK2" => "NTRK2",
    "NTRK3" => "NTRK3",
    "NUTM1" => "NUTM1",
    "P2RY8" => "P2RY8",
    "PALB2" => "PALB2",
    "PARK2" => "PARK2",
    "PARP1" => "PARP1",
    "PARP2" => "PARP2",
    "PARP3" => "PARP3",
    "PAX5" => "PAX5",
    "PBRM1" => "PBRM1",
    "PD-1" => "PDCD1",
    "PDCD1" => "PDCD1",
    "PD-L1" => "CD274",
    "CD274" => "CD274",
    "PD-L2" => "PDCD1LG2",
    "PDCD1LG2" => "PDCD1LG2",
    "PDGFRA" => "PDGFRA",
    "PDGFRB" => "PDGFRB",
    "PDK1" => "PDK1",
    "PIK3C2B" => "PIK3C2B",
    "PIK3C2G" => "PIK3C2G",
    "PIK3CA" => "PIK3CA",
    "PIK3CB" => "PIK3CB",
    "PIK3R1" => "PIK3R1",
    "PIM1" => "PIM1",
    "PMS2" => "PMS2",
    "POLD1" => "POLD1",
    "POLE" => "POLE",
    "PPARG" => "PPARG",
    "PPP2R1A" => "PPP2R1A",
    "PPP2R2A" => "PPP2R2A",
    "PRDM1" => "PRDM1",
    "PRKAR1A" => "PRKAR1A",
    "PRKCI" => "PRKCI",
    "PTCH1" => "PTCH1",
    "PTEN" => "PTEN",
    "PTPN11" => "PTPN11",
    "PTPRO" => "PTPRO",
    "QKI" => "QKI",
    "RAC1" => "RAC1",
    "RAD21" => "RAD21",
    "RAD51" => "RAD51",
    "RAD51B" => "RAD51B",
    "RAD51C" => "RAD51C",
    "RAD51D" => "RAD51D",
    "RAD52" => "RAD52",
    "RAD54L" => "RAD54L",
    "RAF1" => "RAF1",
    "RARA" => "RARA",
    "RB1" => "RB1",
    "RBM10" => "RBM10",
    "REL" => "REL",
    "RET" => "RET",
    "RICTOR" => "RICTOR",
    "RNF43" => "RNF43",
    "ROS1" => "ROS1",
    "RPTOR" => "RPTOR",
    "RSPO2" => "RSPO2",
    "SDC4" => "SDC4",
    "SDHA" => "SDHA",
    "SDHB" => "SDHB",
    "SDHC" => "SDHC",
    "SDHD" => "SDHD",
    "SETD2" => "SETD2",
    "SF3B1" => "SF3B1",
    "SGK1" => "SGK1",
    "SLC34A2" => "SLC34A2",
    "SMAD2" => "SMAD2",
    "SMAD4" => "SMAD4",
    "SMARCA4" => "SMARCA4",
    "SMARCB1" => "SMARCB1",
    "SMO" => "SMO",
    "SNCAIP" => "SNCAIP",
    "SOCS1" => "SOCS1",
    "SOX2" => "SOX2",
    "SOX9" => "SOX9",
    "SPEN" => "SPEN",
    "SPOP" => "SPOP",
    "SRC" => "SRC",
    "STAG2" => "STAG2",
    "STAT3" => "STAT3",
    "STK11" => "STK11",
    "SUFU" => "SUFU",
    "SYK" => "SYK",
    "TBX3" => "TBX3",
    "TEK" => "TEK",
    "TERC" => "TERC",
    "TERT" => "TERT",
    "TET2" => "TET2",
    "TGFBR2" => "TGFBR2",
    "TIPARP" => "TIPARP",
    "TMPRSS2" => "TMPRSS2",
    "TNFAIP3" => "TNFAIP3",
    "TNFRSF14" => "TNFRSF14",
    "TP53" => "TP53",
    "TSC1" => "TSC1",
    "TSC2" => "TSC2",
    "TYRO3" => "TYRO3",
    "U2AF1" => "U2AF1",
    "VEGFA" => "VEGFA",
    "VHL" => "VHL",
    "MMSET" => "NSD2",
    "WHSC1" => "NSD2",
    "NSD2" => "NSD2",
    "WHSC1L1" => "WHSC1L1",
    "WT1" => "WT1",
    "XPO1" => "XPO1",
    "XRCC2" => "XRCC2",
    "ZNF217" => "ZNF217",
    "ZNF703" => "ZNF703");

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
### END OF SUBROUTINES########################################################
############################################################################## 


############
### MAIN ###
############

if(!defined $in or !defined $minPerc) {
    print "\n\tMissing parameter. Try $0 -help to get more information.\n\n";
    exit;
}
    
open(F1, $in) or die "ERROR: Failed to open $in: $!\n";

my %matrix = ();

my $header = <F1>;
my %samples = ();
my %metastasesN = ();
my %block = ();
my @file = ();
while(my $l = <F1>) {
    # split columns 
    my ($id, @a) = split(/\t/, $l);
    my $patient = $id;
    $patient =~ s/(LC\d+\.\d+)\..*/$1/;

    my $block = $id;
    $block =~ s/LC\d+\.\d+\.1\.(\d+)\.1/$1/;
    if(exists $block{$patient}) {
	$block{$patient}->{$id} = keys %{$block{$patient}};
    } else {
	$block{$patient}->{$id} = 0;
    }

    if(exists $metastasesN{$patient}) {
	$metastasesN{$patient}++;
    } else {
	$metastasesN{$patient} = 0;
    }
    push(@file, $l);
}
    
foreach my $l (@file) {
    # split columns 
    my ($id, $snvs, $cnvs, $fusions, $snvsVUS, $cnvsVUS, $fusionsVUS, $tmb, $msStatus, $lohScore, $tumorContentUSZ, $tumorContentFMI) = split(/\t/, $l);

    my $patient = $id;
    $patient =~ s/(LC\d+\.\d+)\..*/$1/;
    
    my @snvs       = split(/,/, $snvs);
    my @cnvs       = split(/,/, $cnvs);
    my @fusions    = split(/,/, $fusions);
    my @snvsVUS    = split(/,/, $snvsVUS);
    my @cnvsVUS    = split(/,/, $cnvsVUS);
    my @fusionsVUS = split(/,/, $fusionsVUS);

    my $t = "";
    if($id =~ /(LC\d+\.\d+)\.1\.1\.1/) {
	$t = "P_";
    } else {
	$t = "M$metastasesN{$patient}"."_".$block{$patient}->{$id}."_";
    }

    $id =~ s/(LC\d+\.\d+)\..*/$1/;
    
    if(defined $vus) {
	push(@snvs, @snvsVUS);
	push(@cnvs, @cnvsVUS);
	push(@fusions, @fusionsVUS);
    }
    # process SNVs. Add them to the matrix hash dependend on their mutation type 
    foreach my $snv (@snvs) {
	if($snv =~ /(.*?) (.*) \(([\d\.]+)%\)/) {
	    my $gene = $fmiGenes{$1};
	    my $mut = $2;
	    my $freq = $3;	    
	    my $mutType = $t."SNV";
	    
	    if(defined $snvSplit) {
		$mutType = $t."Missense";
		$mutType = $t."Nonsense" if($mut =~ /\*/);
		$mutType = $t."Indel" if($mut =~ /ins/ or $mut =~ /del/);
		$mutType = $t."Frameshift" if($mut =~ /fs/);
		$mutType = $t."Splicesite" if($mut =~ /splice/);
		$mutType = $t."Promoter" if($mut =~ /promoter/);
	    }
	    $matrix{$gene}->{$id}->{$mutType} = 1;
	    $samples{$id} = 1;
	} else {
	    print STDERR "WARNING: Unknown SNV format: $id $snv\n";
	}
    }
    # process CNVs
    foreach my $cnv (@cnvs) {
	my $gene = "";
	my $mutType = "";
	if($cnv =~ /(.*) loss/) {
	    $gene = $fmiGenes{$1};
	    $mutType = $t."LOSS";
	} elsif($cnv =~ /(.*?) .*/) {
	    $gene = $fmiGenes{$1};
	    $mutType = $t."AMP";
	} else {
	    print STDERR "WARNING: Unknown CNV format: $id $cnv\n";
	}
	if(!defined $gene){
	    print STDERR "$cnv\n";
	}
	$matrix{$gene}->{$id}->{$mutType} = 1;
	$samples{$id} = 1;
    }
    # process fusions
    foreach my $fusion (@fusions) {
	my $mutType = $t."FUSION";
	if($fusion =~ /(.*)-(.*) fusion/) {
	    my $gene1 = $fmiGenes{$1};
	    my $gene2 = $fmiGenes{$2};

	    if(defined $gene1 and $gene1 ne "N/A") {
		$matrix{$gene1}->{$id}->{$mutType} = 1;
		$samples{$id} = 1;
	    }
	    if(defined $gene2 and $gene2 ne "N/A") {
		$matrix{$gene2}->{$id}->{$mutType} = 1;
		$samples{$id} = 1;
	    }
	} else {
	    print STDERR "WARNING: Unknown Fusion format: $id $fusion\n";
	}
    }
}

# for only taking genes that are mutated in at least $minPerc % of samples
my $nSamples = keys %samples;

# output mutation matrix
# first output header
foreach my $id (sort keys %samples) {
    print "\t$id";
}
print "\n";

foreach my $gene (sort keys %matrix) {
    my $nMutatedSamples = keys %{$matrix{$gene}};

    # only output genes in the minPerc % of samples, except if it is an important melanoma gene
    # /2 as every sample ID are actually two IDs
    if($nMutatedSamples/$nSamples*100/2 >= $minPerc) {

	print $gene;
	foreach my $id (sort keys %samples) {
	    print "\t";
	    if(exists $matrix{$gene}->{$id}) {
		foreach my $mutType (sort keys %{$matrix{$gene}->{$id}}) {
		    print "$mutType;";
		}
	    }
	}
	print "\n";
    }
}
