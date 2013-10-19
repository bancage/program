#!/usr/bin/perl 
#program to change the index of Ref_Seq files
#.fa
@ch=`ls *.fa`;
#print @ch;
for ($chn=0;$chn<=21;$chn++){
#for the Unplaced fa file
if($chn==21){
    if (open INF, "<$ch[$chn]"){
        if (open OUF, ">ch_un.nfa"){
                $c=0;
            while (<INF>){
            	if(/(NW_\d+.\d)/){
                    $c++;
	  	    $ind{"$1"}="NW$c";
		    print OUF ">NW$c","\n";
                }	
   		else{
                    print OUF "$_";
                }
		 
            }
	
         }
        else {
        die ("can not write the  outputfile!");
        }
    }

    else{
        die ("can not open inputfile!");
    }
close INF;
}

#for the placed fa files
else{
    if (open INF, "<$ch[$chn]"){
        @genome=<INF>;
        #print "$genome[$chn]","\n"; 
 	#cat the index
        $genome[0] =~ /(NC_\d+.\d)/;
        #print "$1","\n";
        $NC=$1;
	#print "$NC";
        if($chn==18){
            $chnm="MT";
        }
        elsif($chn==19){
            $chnm="X";
        }
        elsif($chn==20){
            $chnm="Y";
        }
        else{
            $genome[0] =~ /chromosome\s(\d+)/;
	    $chnm=$1;
        }
	#print "$1","\n";
	#write the chr index and it's NC to a hash %ind
	$ind{"$NC"}=$chnm;
	#print %ind;
	$genome[0]=">$chnm\n";
        #print "$genome[$chn]","\n"; 
    }
    else{
        die ("can not open inputfile!");
    }
}

    close INF;
    if(open OUF, ">./chr$chnm.nfa"){
    print OUF join("",@genome);
    }
    else{
        die ("can not write the output file !");
    }
    close OUF; 
}
#print %ind;

#for the gff3 file
    if(open GFF, "<./ref_Sscrofa10.2_scaffolds.gff3"){#ref_Sscrofa10.2_top_level.gff3"){
       if(open NGF, ">./top_level.ngf"){
         while(<GFF>){
         if(/(NW_\d+.\d)/){
             $nnw=$ind{"$1"};
       	     s/NW_\d+.\d/$nnw/;
         }	
         if(/(NC_\d+.\d)/){
             $nnc=$ind{"$1"};
       	     s/NC_\d+.\d/$nnc/;
         
         }	
	     print NGF "$_","\n";
         #else{
         #    print NGF ""
         #}
       }
       }  
      else {
      	  die ("can not write the ngf file!")     
      }
      close NGF;
    
    }

    else{
        die ("can not open the gff3 file !");
    }
    close GFF; 

