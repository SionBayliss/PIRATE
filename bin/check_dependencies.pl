#!/usr/bin/env perl 

# Check dependencies for PIRATE.
# Note: no version check.

# cd-hit
my $cd_hit = 0;
$cd_hit = 1 if `command -v cdhit;`;
$cd_hit = 1 if `command -v cd-hit;`;
print "cd-hit (or alternate invocation cdhit) not found in system path.\n" if $cd_hit == 0;

# blast+
$blast = 0 ;
$blast = 1 if `command -v blastn;`;
print "blast+ not found in system path.\n" if $blast == 0;

# mcl 
my $mcl = 0;
$mcl = 1 if `command -v mcl;`;
print "mcl not found in system path.\n" if $mcl == 0;

# roary (optional)
my $roary = 0;
$roary = 1 if `command -v roary;`;
print "roary not found in system path, cannot use --roary command.\n" if $roary == 0;

# R (optional)
my $R = 0; 
$R = 1 if `command -v R;`;
print "roary not found in system path, cannot use -r command to plot output figures.\n" if $roary == 0;

my $diamond_mkdb = ""; 
my $diamond_bin = "";
$diamond_mkdir = "diamond makedb" if `command -v diamond makedb;`;
$diamond_bin = "diamond blastp" if `command -v diamond blastp;`;

my $diamond_err = 0;
$diamond_err = 1 if $diamond_mkdb eq "";
$diamond_err = 1 if $diamond_bin eq "";
print "Cannot find diamond binaries --diamond cannot be used.\n" if $diamond_err == 1;

# die if dependencies are not available.
if ( ($cd_hit == 0) || ($blast == 0) || ($mcl == 0) ) { die "Dependencies not correctly installed.\n" };

exit
