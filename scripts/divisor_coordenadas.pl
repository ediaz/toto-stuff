#! /usr/bin/perl

#ESTEBAN FERNANDO DÍAZ 12/03/2009

########INTERPOLADOR LINEAL######################


$entrada=$ARGV[0];
$divisor=$ARGV[1];
#$divisor=3;

#BLOQUE DE LECTURA DE PARAMETROS 
print $divisorx;
open(par,"$entrada");
$j=0;
while ( $AR= <par>)
{
$j=$j+1; #j corresponde al horizonte

 	
	$i=$i+1; # corresponde a cada punto del vector de horizontes
	#@values =split /:/, $lee_hor;
	@values = split(';',$AR);
	
	
	
	$x[$j]=@values[0]; #corresponde al valor de la posición en x del punto del horizonte j
	$y[$j]=@values[1]; #corresponde al valor de la posición en z del punto del horizonte j		
	
close FH;

}
$nvalores=$j;
#print $nvalores;


for ($i=1 ; $i<$nvalores ; $i+=1) {
$dx[$i]=($x[$i+1]-$x[$i])/$divisor;
$dy[$i]=($y[$i+1]-$y[$i])/$divisor;
#print "$dx[$i] $x[$i] $x[$i+1]\n" ;
}

#for ($i=$divisor ; $i<($nvalores-1)*$divisor-($divisor-2) ; $i+=1) {
$xxxx=$nvalores*$divisor-($divisor-2);

#for ($i=$divisor ; $i<23 ; $i+=1) {
for ($i=$divisor ; $i<($nvalores)*$divisor+1 ; $i+=1) {
$aux=int($i/$divisor);#+1;
$xx=$aux*$divisor;

$x_int[$i]=$x[$aux]+($i-$xx)*$dx[$aux];

$y_int[$i]=$y[$aux]+($i-$xx)*$dy[$aux];

printf("%.2f\t%.2f\n",$x_int[$i],$y_int[$i] );
}


