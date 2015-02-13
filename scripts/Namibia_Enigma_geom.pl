#! /usr/bin/perl
$pi=3.14159265;


#ESTEBAN FERNANDO DÍAZ 17/04/2009

########INTERPOLADOR LINEAL######################

if ($#ARGV != 0) {
 print "INPUT PARAMETER ERROR!
       
       Usage: 
       ./Namibia_Enigma_geom.pl  [Navigation P190 file] 
       \n";
 exit;
}

$entrada=$ARGV[0];

#$divisor=3;

#BLOQUE DE LECTURA DE PARAMETROS 

open(par,"$entrada");
$j=0;


$stat_sht[0]=659;
$sumX =0;
$sumY =0;
$sumXX =0;
$sumXY =0;
$sumYY =0;
#--------------------------------LECTURA DE SHOTS--------------------------------------------
while ( $AR1= <par>)
{

 	if($AR1=~ m/^(S)/i) {
	$j=$j+1;	

	$FFID[$j]=substr($AR1,21,4);
	if($j==1){
	$FFID[$j-1]=$FFID[$j]-1;
	}
	
	$shot_x[$j]=substr($AR1,47,8); 
	$shot_y[$j]=substr($AR1,55,9);
	$wdepth[$j]=substr($AR1,64,6);
        $sumX =$sumX + $shot_x[$j];
        $sumY =$sumY + $shot_y[$j];
        $sumXX=$sumXX+$shot_x[$j]*$shot_x[$j];
        $sumYY=$sumYY+$shot_y[$j]*$shot_y[$j];
        $sumXY=$sumXY+$shot_x[$j]*$shot_x[$j];
	
		if($j>1){
		$m=$m + ($shot_y[$j]-$shot_y[$j-1])/($shot_x[$j]-$shot_x[$j-1]);
		}
	}	

}

$nshots=$j;
$m=$m/($j-1);

close par;


open(par2,"$entrada");
$station_chan=0;
#$stat_sht[0]=659;
$j=0;
while ( $AR= <par2>)
{
	
 	if($AR=~ m/^(S)/i) {
	$j=$j+1;	
	$stat_sht[$j]=$stat_sht[$j-1]+($FFID[$j]-$FFID[$j-1])*2;
	}	

 	if($AR=~ m/^(R)/i) {
#-------------------------1-----------------------------------------------------------

	$chan=substr($AR,2,3);
	$chan_x=substr($AR,6,8);
	$chan_y=substr($AR,14,9);
	
	$station_x[$stat_sht[$j]-12-$chan]=$station_x[$stat_sht[$j]-12-$chan]+$chan_x;
	$station_y[$stat_sht[$j]-12-$chan]=$station_y[$stat_sht[$j]-12-$chan]+$chan_y;
	$station_fold[$stat_sht[$j]-12-$chan]=$station_fold[$stat_sht[$j]-12-$chan]+1;

#-------------------------2-----------------------------------------------------------
	$chan=substr($AR,28,3);
	$chan_x=substr($AR,32,8);
	$chan_y=substr($AR,40,9);	
	$station_x[$stat_sht[$j]-12-$chan]=$station_x[$stat_sht[$j]-12-$chan]+$chan_x;
	$station_y[$stat_sht[$j]-12-$chan]=$station_y[$stat_sht[$j]-12-$chan]+$chan_y;
	$station_fold[$stat_sht[$j]-12-$chan]=$station_fold[$stat_sht[$j]-12-$chan]+1;

#-------------------------3-----------------------------------------------------------
	$chan=substr($AR,54,3);
	$chan_x=substr($AR,58,8);
	$chan_y=substr($AR,66,9);	
	$station_x[$stat_sht[$j]-12-$chan]=$station_x[$stat_sht[$j]-12-$chan]+$chan_x;
	$station_y[$stat_sht[$j]-12-$chan]=$station_y[$stat_sht[$j]-12-$chan]+$chan_y;
	$station_fold[$stat_sht[$j]-12-$chan]=$station_fold[$stat_sht[$j]-12-$chan]+1;
	}	

}

$station_fold[$stat_sht[$j]-12-$chan];
close par2;
for ($i=1 ; $i<$#station_fold +1 ; $i+=1) {

	$station_x[$i]=$station_x[$i]/$station_fold[$i];
	$station_y[$i]=$station_y[$i]/$station_fold[$i];
}



###################################TERMINA BLOQUE DE COORDENADAS REALES ##############################################

##################################EMPIEZAN LAS COORDENADAS FICTICIA ##################################################





$a2=($shot_y[$nshots]-$shot_y[1])/($shot_x[$nshots]-$shot_x[1]);
if($shot_x[$nshots]>$shot_x[1]){
$crece=1;
}
else{
$crece=0;
}

$a=$m;

print "pendiente=  $a\n";
#--------------------------------STATION CALULADO-----------------------------
if($crece==1){
$station_x_calc[1]=$shot_x[1]-(165+647*12.5)*cos(atan2($a,1));
$station_y_calc[1]=$shot_y[1]-(165+647*12.5)*sin(atan2($a,1));
}
else{
$station_x_calc[1]=$shot_x[1]+(165+647*12.5)*cos(atan2($a,1));
$station_y_calc[1]=$shot_y[1]+(165+647*12.5)*sin(atan2($a,1));
}	
$rec=join "","rec_",$entrada,".txt";

open(out_stat,">$rec");
$i=1;
printf out_stat ("%5s%10.1f%10.1f%10.1f%10.1f\n",$i,$station_x_calc[$i],$station_y_calc[$i],$station_x[$i],$station_y[$i]);
for ($i=2 ; $i<2*$nshots +661 ; $i+=1) {
	if($crece==1){
	$station_x_calc[$i]=$station_x_calc[$i-1]+12.5*cos(atan2($a,1));
	$station_y_calc[$i]=$station_y_calc[$i-1]+12.5*sin(atan2($a,1));
	}
	else{
	$station_x_calc[$i]=$station_x_calc[$i-1]-12.5*cos(atan2($a,1));
	$station_y_calc[$i]=$station_y_calc[$i-1]-12.5*sin(atan2($a,1));
	}
printf out_stat ("%5s%10.1f%10.1f%10.1f%10.1f\n",$i,$station_x_calc[$i],$station_y_calc[$i],$station_x[$i],$station_y[$i]) ;
}





close out_stat ;

##################################################################################################

#-------------------------------SHOT CALCULADO------------------------------------
$sht=join "","shot_",$entrada,".txt";
open(shot_out,">$sht");

$shot_x_calc[1]=$shot_x[1];
$shot_y_calc[1]=$shot_y[1];
$i=1;

printf shot_out ("%5s%5s%5s%10.1f%10.1f%10.1f%10.1f%10.1f\n",$i,$FFID[$i],$stat_sht[$i],$shot_x_calc[$i],$shot_y_calc[$i],$shot_x[$i],$shot_y[$i],$wdepth[$i]) ;
for ($i=2 ; $i<$#stat_sht +1 ; $i+=1) {
	if($crece==1){
	$shot_x_calc[$i]=$shot_x_calc[$i-1]+25.0*cos(atan2($a,1));
	$shot_y_calc[$i]=$shot_y_calc[$i-1]+25.0*sin(atan2($a,1));
	}
	else{
	$shot_x_calc[$i]=$shot_x_calc[$i-1]-25.0*cos(atan2($a,1));
	$shot_y_calc[$i]=$shot_y_calc[$i-1]-25.0*sin(atan2($a,1));
	}
printf shot_out ("%5s%5s%5s%10.1f%10.1f%10.1f%10.1f%10.1f\n",$i,$FFID[$i],$stat_sht[$i],$shot_x_calc[$i],$shot_y_calc[$i],$shot_x[$i],$shot_y[$i],$wdepth[$i]) ;
}
	
close shot_out;




#############################################RELATION ###############################################

$rel=join "","rel_",$entrada,".txt";
open(rel_out,">$rel");

for ($i=1 ; $i<$#stat_sht +1 ; $i+=1) {
printf rel_out ("%5s%5s%5s%5s%5s%5s\n",$i,$stat_sht[$i],1,$stat_sht[$i]-13,648,$stat_sht[$i]-660) ;
}
close rel_out;
