#! /usr/bin/perl


if ($#ARGV != 4) {
 print "INPUT PARAMETER ERROR!
       
#DEVELOPED BY : ESTEBAN DÍAZ PANTIN 17/04/2009
#
#Usage
#
#Type in the command window the following: 
#
#./LEPLAC_geom.pl [navigation p190 file] [number of channels] [minimum offset] [station distance] [shot distance]
#
#where nav.p190 is the P1/90 file to be process
#
#
#
#
#This program does a geometry for 2D marine data, from P190 navigation
#data.
#
#The program prints: receiver table, shot table and relation table.
#
#The receiver table has:
#Station X_calc Y_calc mean_X mean_Y
#
#The shot table has:
#Shot(seq) FFID shot_stat(the nearest) X_calc Y_calc X_original Y_original Wdepth
#
#The relation table has:
#Shot(seq) shot_stat(the nearest) chan1 station_chan1 last_chan station_last_chan
#
#
#IMPORTANT: this script is configurated to read p190 STANDARD data, normally should be OK
#, but in some cases the user has to change the position and offset for each variable to 
# be read. This variables are: 
#From shot: shot_x, shot_y, FFID, wdepth
#From rec: chan1x, chan1y, chan2x, chan2y, chan3x, chan3y (the p190 has three collumns of
#receiver coordinates data.
";
 exit;
}





$entrada=$ARGV[0];
$pi=3.14159265;
#$divisor=3;

#BLOQUE DE LECTURA DE PARAMETROS 

open(par,"$entrada");
$j=0;


#====================================================================
#INPUT PARAMETERS:
$ncanales=$ARGV[1]; #408; #NUMBER OF CHANNELS OF CABLE

$offset_minimo=$ARGV[2]; #116; #MINIMUM OFFSET

$delta_stat=$ARGV[3]; #25; #DISTANCE BETWEEN REC STATIONS

$delta_shot=$ARGV[4]; #50; #DISTANCE BETWEEN SHOT STATIONS

#====================================================================





$stat_sht[0]=$ncanales+int($offset_minimo/$delta_stat)-int($delta_shot/$delta_stat);

$shot_station_step=int($delta_shot/$delta_stat);


$sumX =0;
$sumY =0;
$sumXX =0;
$sumXY =0;
$sumYY =0;
#--------------------------------LECTURA DE SHOTS--------------------------------------------
while ( $AR1= <par>)
{
 #j corresponde al horizonte

	
	
 	if($AR1=~ m/^(S)/i) {
	$j=$j+1;	

	$FFID[$j]=substr($AR1,20,5);
           	if($j==1){
	     	$FFID[$j-1]=$FFID[$j]-1;
	   	}
	
	$shot_x[$j]=substr($AR1,46,9); 
	$shot_y[$j]=substr($AR1,55,9);
	$wdepth[$j]=substr($AR1,64,6);
#print "$wdepth[$j]\n";
        $sumX =$sumX + $shot_x[$j];
        $sumY =$sumY + $shot_y[$j];
        $sumXX=$sumXX+$shot_x[$j]*$shot_x[$j];
        $sumYY=$sumYY+$shot_y[$j]*$shot_y[$j];
        $sumXY=$sumXY+$shot_y[$j]*$shot_x[$j];
	
		if($j>1){
		$m=$m + ($shot_y[$j]-$shot_y[$j-1])/($shot_x[$j]-$shot_x[$j-1]);
		}
	#printf("%10s%15.3f%15.3f\n",$FFID[$j],$shot_x[$j],$shot_y[$j]);	
	}	

}

$nshots=$j;

#LINEAR REGRETION FOR INCLINATION CALCULATION:
$m_reg=($sumX*$sumY - $nshots*$sumXY)/($sumX*$sumX - $nshots*$sumXX); 

$m=$m/($j-1);





open(par2,"$entrada");
$station_chan=0;

$j=0;
while ( $AR= <par2>)
{
 #j corresponde al horizonte

	
	
 	if($AR=~ m/^(S)/i) {
	$j=$j+1;	
	
	
	$stat_sht[$j]=$stat_sht[$j-1]+int($delta_shot/$delta_stat);
#	print "$stat_sht[$j] $j\n";
	}	

 	if($AR=~ m/^(R)/i) {
#-------------------------1-----------------------------------------------------------

	$chan=substr($AR,2,3);
	$chan_x=substr($AR,5,9);
	$chan_y=substr($AR,14,9);
        $channel_station=$stat_sht[$j]-(int($offset_minimo/$delta_stat)-1)-$chan;
#print "CHANNEL STAT: $channel_station\n";
#print "1: $chan_x  $chan_y\n";
	$station_x[$channel_station]=$station_x[$channel_station]+$chan_x;
	$station_y[$channel_station]=$station_y[$channel_station]+$chan_y;
	$station_fold[$channel_station]=$station_fold[$channel_station]+1;

#-------------------------2-----------------------------------------------------------
	$chan=substr($AR,28,3);
	$chan_x=substr($AR,31,9);
	$chan_y=substr($AR,40,9);	
        $channel_station=$stat_sht[$j]-(int($offset_minimo/$delta_stat)-1)-$chan;

#print "2: $chan_x  $chan_y\n";
#print "CHANNEL STAT: $channel_station\n";
	$station_x[$channel_station]=$station_x[$channel_station]+$chan_x;
	$station_y[$channel_station]=$station_y[$channel_station]+$chan_y;
	$station_fold[$channel_station]=$station_fold[$channel_station]+1;
#-------------------------3-----------------------------------------------------------
	$chan=substr($AR,54,3);
	$chan_x=substr($AR,57,9);
	$chan_y=substr($AR,66,9);	
        $channel_station=$stat_sht[$j]-(int($offset_minimo/$delta_stat)-1)-$chan;
#print "3: $chan_x  $chan_y\n";
#print "CHANNEL STAT: $channel_station\n";
	$station_x[$channel_station]=$station_x[$channel_station]+$chan_x;
	$station_y[$channel_station]=$station_y[$channel_station]+$chan_y;
	$station_fold[$channel_station]=$station_fold[$channel_station]+1;
	}	

}


close par2;
#open(out_stat,">$rec");
for ($i=1 ; $i<$#station_fold +1 ; $i+=1) {

#print "$i $station_fold[$i] \n";
	$station_x[$i]=$station_x[$i]/$station_fold[$i];
	$station_y[$i]=$station_y[$i]/$station_fold[$i];
#printf ("%10s %10.1f %10.1f\n",$i,$station_x[$i],$station_y[$i]) ;
}
#close out_rec;

close par;



####################TERMINA BLOQUE DE COORDENADAS REALES ##################################

#####################EMPIEZAN LAS COORDENADAS CALCULADAS###################################





$a2=($shot_y[$nshots]-$shot_y[1])/($shot_x[$nshots]-$shot_x[1]);
#print "a: $m   a2:  $a2  \n";


$a=$m_reg;

print "pendiente=  $a\n";
#--------------------------------STATION CALULADO-----------------------------


#Checking the direction of vessel adquisition: 
#NE
if($shot_x[2]>$shot_x[1] && $shot_y[2]>$shot_y[1]){$dir_x=1;$dir_y=1; };
#SW
if($shot_x[2]<$shot_x[1] && $shot_y[2]<$shot_y[1]){$dir_x=-1;$dir_y=-1; };
#NW
if($shot_x[2]<$shot_x[1] && $shot_y[2]>$shot_y[1]){$dir_x=-1;$dir_y=1; };
#SE
if($shot_x[2]>$shot_x[1] && $shot_y[2]<$shot_y[1]){$dir_x=1;$dir_y=-1; };





$station_x_calc[1]=$shot_x[1]-($offset_minimo+($ncanales-1)*$delta_stat)*$dir_x*abs(cos(atan2($a,1)));
$station_y_calc[1]=$shot_y[1]-($offset_minimo+($ncanales-1)*$delta_stat)*$dir_y*abs(sin(atan2($a,1)));


#WRITTING RECEIVER TABLE

$rec=join "","rec_",$entrada,".txt";

open(out_stat,">$rec");
$i=1;
printf out_stat ("%5s%10.1f%10.1f%10.1f%10.1f\n",$i,$station_x_calc[$i],$station_y_calc[$i],$station_x[$i],$station_y[$i]);
for ($i=2 ; $i<int($delta_shot/$delta_stat)*$nshots + ($stat_sht[0]+2)+ int($offset_minimo/$delta_stat)-1; $i+=1) {

	$station_x_calc[$i]=$station_x_calc[$i-1]+$dir_x*$delta_stat*abs(cos(atan2($a,1)));
	$station_y_calc[$i]=$station_y_calc[$i-1]+$dir_y*$delta_stat*abs(sin(atan2($a,1)));
 if($station_x[$i]==0){
 $station_x[$i]=$station_x[$i-1]+$dir_x*$delta_stat*abs(cos(atan2($a,1)));
 $station_y[$i]=$station_y[$i-1]+$dir_y*$delta_stat*abs(sin(atan2($a,1))); 
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
#WRITTING SHOT TABLE
printf shot_out ("%5s%8s%5s%10.1f%10.1f%10.1f%10.1f%10.0f\n",$i,$FFID[$i],$stat_sht[$i],$shot_x_calc[$i],$shot_y_calc[$i],$shot_x[$i],$shot_y[$i],$wdepth[$i]) ;
for ($i=2 ; $i<$#stat_sht +1 ; $i+=1) {

	$shot_x_calc[$i]=$shot_x_calc[$i-1]+$dir_x*$delta_shot*abs(cos(atan2($a,1)));
	$shot_y_calc[$i]=$shot_y_calc[$i-1]+$dir_y*$delta_shot*abs(sin(atan2($a,1)));

printf shot_out ("%5s%8s%5s%10.1f%10.1f%10.1f%10.1f%10.1f\n",$i,$FFID[$i],$stat_sht[$i],$shot_x_calc[$i],$shot_y_calc[$i],$shot_x[$i],$shot_y[$i],$wdepth[$i]) ;
}
	
close shot_out;


#WRITTING RELATION TABLE

#############################################RELATION ###############################################

$rel=join "","rel_",$entrada,".txt";
open(rel_out,">$rel");

for ($i=1 ; $i<$#stat_sht +1 ; $i+=1) {
printf rel_out ("%5s%5s%5s%5s%5s%5s\n",$i,$stat_sht[$i],1,$stat_sht[$i]-int($offset_minimo/$delta_stat),$ncanales,$stat_sht[$i]-($ncanales+int($offset_minimo/$delta_stat))+1 ) ;
}
close rel_out;













