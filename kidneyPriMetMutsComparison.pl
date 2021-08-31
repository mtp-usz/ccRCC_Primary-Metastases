#!/usr/bin/perl -w
# Author: Abdullah Kahraman
# Date: 26.08.2021

###############################################################################
###############################################################################
### Lists mutations from KidneyPriMet project line by line.                 ###
###############################################################################
###############################################################################

use strict;
use warnings;
use Getopt::Long;
no autovivification;

my (
    # variable for parameters which are read in from commandline
    $help,
    $vus,
    $all,
    );

##############################################################################
### read all needed parameters from commandline ##############################

&GetOptions(
    "help!"      => \$help,    # print this help
    "vus!"       => \$vus,     # print out VUS only
    "all!"       => \$all,     # print out all mutations
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

###############################################################################
sub analyse {
    my ($h, $p, $title) = @_;

    foreach my $lc (sort keys %$h) {    
	for(my $i = 2; $i < 10; $i++) {
	    if(exists $h->{$lc}->{1} and exists $h->{$lc}->{$i} and
	       exists $p->{$lc}->{1} and exists $p->{$lc}->{$i} and $p->{$lc}->{1} ne $p->{$lc}->{$i}) {
		my @pri = split(/,/, $h->{$lc}->{1});
		my @met = split(/,/, $h->{$lc}->{$i});

		next if(@pri == 0 and @met == 0);
		if(@pri > 0) {
		    foreach my $mut (@pri) {
			my $gene = $mut;
			$gene =~ s/ .*//;
			print "$lc\tPrimary\t$title\t$gene\t$mut\n";
		    }
		} else {
			print "$lc\tPrimary\t$title\t-\t-\n";
		}
		if(@met > 0) {
		    foreach my $mut (@met) {
			my $gene = $mut;
			$gene =~ s/ .*//;
			print "$lc\tMetastase-$i\t$title\t$gene\t$mut\n";
		    }
		} else {
			print "$lc\tMetastase-$i\t$title\t-\t-\n";
		}
	    }
	}
    }
}
##############################################################################
### END OF SUBROUTINES########################################################
############################################################################## 


############
### MAIN ###
############

my (%s, %c, %r, %p);
while(my $l = <>) {
    chomp($l);
    my @a = split(/\t/, $l);
    next if($a[3] !~ /^USZ/);
    my $id = $a[0];
    my $sv = $a[10];
    my $cna = $a[11];
    my $rearr = $a[12];
    my $svVUS = $a[13];
    my $cnaVUS = $a[14];
    my $rearrVUS = $a[15];
    my $specimen = $a[8];

    if(defined $vus) {
	$sv = $svVUS;
	$cna = $cnaVUS;
	$rearr = $rearrVUS;
    } elsif(defined $all){
	$sv .= ",$svVUS";
	$cna .= ",cnaVUS";
	$rearr .= ",$rearrVUS";
    }

    next if(
	$id =~ /LC2019\.224\./ or
	$id =~ /LC2019\.259\./ or
	$id =~ /LC2019\.261\./ or
	$id =~ /LC2019\.268\./ or
	$id =~ /LC2019\.269\./ or
	$id =~ /LC2019\.282\./ or
	$id =~ /LC2019\.288\./ or		
	$id =~ /LC2019\.289\./ or
	$id =~ /LC2020\.504\./);
    my $lc = $id;
    $lc =~ s/(LC20\d\d\.\d+)\..*/$1/;
    my $block = $id;
    $block =~ s/LC20\d\d\.\d+\.\d+\.(\d+)\.\d+/$1/;

    $s{$lc}->{$block} = $sv;
    $c{$lc}->{$block} = $cna;
    $r{$lc}->{$block} = $rearr;
    $p{$lc}->{$block} = $specimen;
}

    &analyse(\%s, \%p, "SV");
    &analyse(\%c, \%p, "CNA");
    &analyse(\%r, \%p, "Fusion");

