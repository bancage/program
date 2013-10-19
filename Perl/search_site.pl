#!/usr/bin/perl -w
#program to search the binding site of our interested transcription factor from TRANSFAC 
#output file "site.dat", by JessicaBan Jun 27th,2013.use strict;
use strict;
my $infile = "./site.dat";
my $outfile = "./chickhif1a.txt";
my $TF = "(HIF-?1 ?alpha)";
my $spe = "chick";
local $/ = "//";
open my $inf, "<$infile" or die "$infile is no exists";
open my $ouf, ">$outfile" or die "$outfile can not be written";
while(<$inf>){
    if (m/^OS.*$spe/im && m/.*$TF/im){
           print $ouf $_;
       }
}
