#!/usr/bin/perl -w
#program to count the HRE 
$nseq = "./n_hre.txt";
$e_name = "HRE";
$e_seq1 = "ACGTG";    #"RCGTG"
$e_seq2 = "GCGTG";    #"RCGTG"
$chn_fr = 1;
$chn_to = 21;
for ($fn=$chn_fr;$fn<=$chn_to;$fn++){


#$infile = "./test$fn.txt";
$infile = "./ssc_ref_Sscrofa10_chr$fn.fa.bak";

if (open INF, "<$infile"){
    @genome = <INF>;
}
else {
    die ("can't open genome file\n");
}

#close the input file
close  INF;    

#join them
$gn = join('',@genome);    
chomp($gn);
$gn =~ s/\s+//g;

#count the number of the element
#$n_e[$fn]=$gn =~ /($e_seq1)|($e_seq2)/g;
while($gn =~ /($e_seq1)|($e_seq2)/g){$n_e[$fn]++;}

#write the number of the element into the outfile
if (open NS, ">>$nseq"){
    print NS "element number of Chr $fn is: $n_e[$fn]", "\n";
    if ($fn==$chn_to){
        $s = &sum(@n_e);
        print NS "The genomic element number is: $s !", "\n" ;
    }
}
else {
    die ("can't append the number file\n");
}

#close the input file
close  NS;    
}

sub sum{
    foreach (@_){
        $su += $_;
    }    
        $su;
}
