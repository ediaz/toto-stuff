#!/usr/bin/perl 
#
#**********************************************************
#Este script selecciona los shots dentro de un 
#area rectangular con rumbo oblicuo a partir de 
#archivos de navegación UKOOA90
#sintanxis: 
#extrate_area.pl <archivo de entrada> <archivo de salida>
#  

#*************************************************************
#     / \
#    / a \
#   /     \
#  /       \
# 3\b       \  
#   \        \
#    \      b/2
#     \     /
#      \   /
#       \a/
#        1 
# a=b=90 degrees
#************************************************************

#El programa rota el sistema de coordenadas de forma tal que
#posterior a la rotacion quedaria de la siguiente forma:
# 
#            dx
#	  ________
#	 |b      a|
#	 |        |
#	 |        |
#	 |        |dy
#	 |        |
#	 |        |
#	 |a______b|
#       /\    
#      /||\
#       || 
#     (xo,yo)
#
# Lo que hace mas facil la selección de los puntos

#dx= d 1==>2
#dy= d 1==>3


#punto (0,0):
$xo=584518.903;
$yo=4089374.88;
#--------------
$pi=3.14159265;
$dx=5000; 
$dy=6000;
$iii=0;

#theta: ángulo de rumbo del área (medido con respecto al Norte y con sentido horario)
$theta=39.7012;
$linea_S=0;
$gap_full_fold=1000;
$entrada="$ARGV[0]";
$limite_y1=0;
$limite_y2=0;
open (FILE, $entrada);
$salida="ARGV[1]";
open (out,">$salida");
$n_linea=0;
while ($linea=<FILE>)
{

$dir=1;

$n_linea=$n_linea+1;

	if($linea =~ m/^(H)/i) {
	print out $linea;
	}
$shot="no escribir";
	if($linea =~ m/^(S)/i) {


		$sail_line=substr($linea,6,7);
		$A=substr($linea,0,20);
		$B=substr($linea,21,79);
		$x=substr($linea,47,8);
		$y=substr($linea,55,9);

		$x_rotado=($x-$xo)*cos($theta*$pi/180)+($y-$yo)*sin($theta*$pi/180);
		$y_rotado=($y-$yo)*cos($theta*$pi/180)-($x-$xo)*sin($theta*$pi/180);


			if($sail_line=="946P030" || $sail_line=="934P032" || $sail_line=="922P034" || $sail_line=="922F036" || $sail_line=="922F038" ) {
			$dir=-1;
			}
# dir es el sentido de la linea de navegación
			if($dir == 1) {
			$limite_y2=6050+$gap_full_fold+120;
			$limite_y1=-50;
			}

#dir = 1 denota dirección increase
#dir =-1 denota dirección decrease

			if($dir == -1) 
			{
			$limite_y1=-50-$gap_full_fold-120;
			$limite_y2=6050;
			}
#print "$dir $limite_y1 $limite_y2\n";
			if($x_rotado<5050 && $x_rotado>-50.0 && $y_rotado< $limite_y2 && $y_rotado>$limite_y1){
			$linea_S=$n_linea;
			$iii=$iii+1;
			$shot="escribir";
			}
			else
			{
			$shot="no escibir"
			}
	}

	if($linea_S==$n_linea){
	print out "$A $B";
	}
	else
	{
		if($linea_S>0 && ($n_linea-$linea_S)<(289)){
		print out $linea;
		}
	}
}
print "Puntos: $iii\n";
close out;
