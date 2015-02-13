#! /usr/bin/perl 

$entrada=$ARGV[0];
$busca=int($ARGV[1]);
$separador=$ARGV[2];

open(par,"$entrada");
$j=0;
$line_t=999999;


$shot_hdr="SHOT INFORMATION:                                                               
S     PT        X          Y         DS      UH    ELEV    DATUM    LGWS   ELCOR
     SHWS      TOTAL      VREP       D1      V1     D2      V2       D3      V3 
      D4        V4         D5        V5     VSUB   FOLD        
";

print $shot_hdr ;
               
while ($AR= <par>)
{

$line=$line+1;
$X=$AR;
$X=~ s/S//g;
@values = split(/\ +/,$X);

if($AR =~ m/^(S)/i && (@values[1]!="PT")){
if( @values[1]<(($busca+2)*$separador) && @values[1]=>(($busca+1)*$separador)) {
$corr=@values[1]-(($busca+1)*$separador);
#print "@values[1] $corr\n";
$X=~ s/@values[1]/     /g;
#print "$X";
$xx=substr($AR,16,80);
	if($corr< $separador && $corr>0 ){
		printf ("S%10i     %s",$corr,$xx) ;
	$line_t=$line;	
	}
}




}
#$j=$j+1;

if($line==$line_t+1 ||$line==$line_t+2){print $AR}



#}
}
close par;

#########################################################REC
$rec_hdr="END OF SHOT INFORMATION.                                                        
RECEIVER INFORMATION:                                                           
R     PT        X          Y                       ELEV    DATUM    LGWS   ELCOR
     SHWS      TOTAL      VREP       D1      V1     D2      V2       D3      V3 
      D4        V4         D5        V5     VSUB   FOLD                         
";

open(par,"$entrada");
$j=0;
$line_t=999999;



print $rec_hdr ;
               
while ($AR= <par>)
{

$line=$line+1;
$X=$AR;
$X=~ s/R//g;
@values = split(/\ +/,$X);

if($AR =~ m/^(R)/i && (@values[1]!="PT")){
if( @values[1]<(($busca+2)*$separador) && @values[1]=>(($busca+1)*$separador)) {
$corr=@values[1]-(($busca+1)*$separador);
#print "@values[1] $corr\n";
$X=~ s/@values[1]/     /g;
#print "$X";
$xx=substr($AR,14,80);
	if($corr< $separador && $corr>0 ){
		printf ("R%10i   %s",$corr,$xx) ;
	$line_t=$line;	
	}
}




}
#$j=$j+1;

if($line==$line_t+1 ||$line==$line_t+2){print $AR}



#}
}
print "END OF RECEIVER INFORMATION.
";
close par;
