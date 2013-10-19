#!/usr/bin/perl -w
#program to do the enzyme digestion 
$nseq = "./n_seq.txt";
$re_name = "MspI";
$re_seq = "CCGG";
$chn_fr = 1;
$chn_to = 21;
for ($fn=$chn_fr;$fn<=$chn_to;$fn++){

$infile = "./ssc_ref_Sscrofa10_chr$fn.fa";
$outfile = "./ssc_ref_Sscrofa10_chr$fn.out.fa";
#read the input file into @genome
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
#print $gn, "\n";

#cut them
@list = split ($re_seq, $gn);    

#calculate the number of cutted sequences and write the number into the n_seq.txt file
$n_seq[$fn] = $#list + 1;
if (open NS, ">>$nseq"){
    print NS "cutted seq number of Chr $fn is: $n_seq[$fn]", "\n";
    if ($fn==$chn_to){
        $s = &sum(@n_seq);
        print NS "The genomic cutted seq number is: $s !", "\n" ;
    }
}
else {
    die ("can't append the number file\n");
}

#close the input file
close  NS;    

#make them up
#the first sequence
$list[0] .= "C";

#the last sequence
$list[$n_seq-1] = "CGG".$list[$n_seq-1];

#the middle ones
for ($i = 1; $i <$#list; $i++){
    $list[$i] = "CGG".$list[$i]."C"; 
    # print $i;
    #print $list[$i], "\n";
    }
#output the resulting sequnces into *_out.fa
if (open OUF, ">$outfile"){
        foreach $list(@list){
            print OUF $list,"\n";
        }
     
}
else {
    die ("can't creat output file\n");
}

#close the output file
close OUF;   
}


sub sum{
    foreach (@_){
        $su += $_;
    }    
        $su;
}
