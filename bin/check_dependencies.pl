#!/usr/bin/env perl 

use strict;
use warnings;

# Check dependencies for PIRATE.
# Note: currently no version check.

# parallel
my $parallel = 0;
$parallel = 1 if `command -v parallel;`;
print " - ERROR: GNU parallel not found in system path.\n" if $parallel == 0;

# cd-hit
my $cd_hit = 0;
$cd_hit = 1 if `command -v cdhit;`;
$cd_hit = 1 if `command -v cd-hit;`;
print " - ERROR: cd-hit (or alternate invocation cdhit) not found in system path.\n" if $cd_hit == 0;

# blast+
my $blast = 0 ;
$blast = 1 if `command -v blastn;`;
print " - ERROR: blast+ not found in system path.\n" if $blast == 0;

# mcl 
my $mcl = 0;
$mcl = 1 if `command -v mcl;`;
print " - ERROR: mcl not found in system path.\n" if $mcl == 0;

# die if dependencies are not available.
if ( ($cd_hit == 0) || ($blast == 0) || ($mcl == 0) || ($parallel == 0) ) { die "\n - ERROR: Dependencies not correctly installed.\n" };

# check optional dependencies

# diamond
my $diamond_mkdb = ""; 
my $diamond_bin = "";
$diamond_mkdb = "diamond makedb" if `command -v diamond makedb;`;
$diamond_bin = "diamond blastp" if `command -v diamond blastp;`;

my $diamond_err = 0;
$diamond_err = 1 if $diamond_mkdb eq "";
$diamond_err = 1 if $diamond_bin eq "";
print "\n - WARNING: cannot find diamond binaries, cannot use --diamond command.\n" if $diamond_err == 1;

# R
my $R = 0; 
$R = 1 if `command -v R;`;
print " - WARNING: R not found in system path, cannot use -r command.\n" if $R == 0;

# fasttree
my $ft = 0; 
$ft = 1 if `command -v fasttree;`;
print " - WARNING: fasttree not found in system path, a binary presence-absence tree will not be created.\n" if $ft == 0;

exit
