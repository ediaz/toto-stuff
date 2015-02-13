#! /usr/bin/perl
$pi=3.14159265;


#ESTEBAN FERNANDO DÍAZ 17/04/2009

$entrada=$ARGV[0];
$entrada1=$ARGV[1];
#$divisor=3;

#BLOQUE DE LECTURA DE PARAMETROS 

open(par,"$entrada");
$j=0;

$ncanales=408;

$offset_minimo=116;

$delta_stat=25;

$delta_shot=50;

$stat_sht[0]=$ncanales+int($offset_minimo/$delta_stat)-int($delta_shot/$delta_stat);

$shot_station_step=2;


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

	$FFID[$j]=substr($AR1,19,4);
	if($j==1){
	$FFID[$j-1]=$FFID[$j]-1;
	}
	
	$shot_x[$j]=substr($AR1,46,10); 
	$shot_y[$j]=substr($AR1,55,10);
	
#print "$wdepth[$j]\n";
        $sumX =$sumX + $shot_x[$j];
        $sumY =$sumY + $shot_y[$j];
        $sumXX=$sumXX+$shot_x[$j]*$shot_x[$j];
        $sumYY=$sumYY+$shot_y[$j]*$shot_y[$j];
        $sumXY=$sumXY+$shot_y[$j]*$shot_x[$j];
	
		
	#printf("%10s%15.3f%15.3f\n",$FFID[$j],$shot_x[$j],$shot_y[$j]);	
	}	

}

$nshots=$j;
$m_reg=($sumX*$sumY - $nshots*$sumXY)/($sumX*$sumX - $nshots*$sumXX);




$j=0;

open(par2,"$entrada1");

while ( $AR1= <par2>)
{
 #j corresponde al horizonte

	
	
 	if($AR1=~ m/^(R)/i) {
	$j=$j+1;	

	$stat[$j]=substr($AR1,19,4);
	if($j==1){
	$FFID[$j-1]=$FFID[$j]-1;
	}
	
	$stat_x[$j]=substr($AR1,46,10); 
	$stat_y[$j]=substr($AR1,55,10);
		
	
	}	

}

$nstats=$j;



for ($i=1; $i<$nshots+1;$i++){
#loop of source	
	for ($j=1; $j<$nstats+1;$j++){
		#loop of station
	$cmpx[$i][$j]=($shot_x[$i]+$stat_x[$j])/2;
	$cmpy[$i][$j]=($shot_y[$i]+$stat_y[$j])/2;	
#print "$cmpx[$i][$j] $cmpy[$i][$j] \n";
	}

}

$ncpd=$nstats+$nshots-1;




$k=1;



for ($k=1; $k<445 ; $k++){
$i=$k;
$j=1;
$f=0;


 #   while ($i+$j==$k+1 && $i>0 && $j>0){
 print "=================== $k===============\n";
    while ( $i>0 && $j>0){
   
    $cmp_x[$k] = $cmp_x[$k]+($shot_x[$i]+$stat_x[$j])/2;
    $cmp_y[$k] = $cmp_y[$k]+($shot_y[$i]+$stat_y[$j])/2;
    $f=$f+1;
	$x=($shot_x[$i]+$stat_x[$j])/2;
	$y=($shot_y[$i]+$stat_y[$j])/2;
    print "$i $j  $x $y\n";
    $i=$i-1;
    $j=$j+1;
     }

$cmp_x[$k]=$cmp_x[$k]/$f;
$cmp_y[$k]=$cmp_y[$k]/$f;
$fold[$k]=$f;
#printf("%8d%9.1f%9.1f%8d\n",$k,$cmp_x[$k],$cmp_y[$k],$fold[$k]);

}

