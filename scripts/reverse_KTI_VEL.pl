#!/usr/bin/perl 

#Reverse the KTI vel order of CDP based on the first 
#and last CDP

#USAGE MESSAGE
if ($#ARGV != 1) {
 print "INPUT PARAMETER ERROR!
       
       Usage: 
       ./reverse_KTI_VEL.pl  [vel or mute file] [cdp_min+cdp_max] > [output file]
       \n";
 exit;
}




$entrada="$ARGV[0]";
$n="$ARGV[1]";

open (FILE, "<$entrada");

while (<FILE>){

@cdp_split= split(/\ +/,$_);
$test=true;
   if( /CDP/ || /LOC/ ) #Searching for CDP (vel files) or LOC (mute files)
    { 
     $cdp_in=@cdp_split[2];
  
     $cdp_out=$n-$cdp_in;
     
     $print_out=join "",$cdp_out,"   < Old value = ",$cdp_in," >";    

     s/$cdp_in/$print_out/g;
    }
      print $_;
    	
}
