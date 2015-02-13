#!/usr/bin/perl 

#FORMATO VEL FOCUS A FORMATO VEL ANP
#En la linea de comando escribir:
#focus2anp.pl [nombre de archivo de velocidades focus] [nombre de la linea]
#el nombre de la linea no puede tener mas de 11 caracteres
#para salvarlo en un archivo escribir:
#focus2anp.pl [nombre del archivo de velocidades focus] [nombre de la linea] > [nombre de archivo de salida]


$entrada="$ARGV[0]";
$name="$ARGV[1]";

open (FILE, $entrada);

$n_linea=0;
$sample=1;
$cdp_ant=-9999999;
$cdp=-99999999;
$h=0;
$elementos=0;
print "LINE       $name\n";
$imprimibles=0;
while ($linea=<FILE>)
{
$n_linea=$n_linea+1;


   if($linea =~ m/^H/i) #lineas que comienzan con handvel
    {
    $h=$h+1;
    @cdp_split= split(/\ +/,$linea);
    $cdp_ant=$cdp;
    
    $cdp=@cdp_split[1];
    if($cdp>0){
    $imprimibles=$imprimibles+1;
    }
#   if($h==1){
#   printf("SPNT%22s\n",$cdp);
#    }
    if($imprimibles==1 && $cdp_ant<0 && $cdp>0){
    printf("SPNT%21s\n",$cdp);
    }
  	
    $i=0;
    $n=0;
             
         
     if($h<=>1 && $cdp>0 && $cdp_ant>0 ){
#Dentro de este if ya tengo todos los elementos de la funcion de velocidad que corresponde
#al cdp $cdp


     printf ("VELF%16s");
     
          for ($i=0;$i<$elementos;$i++) {
              printf ("%5d",@tvpairs[$i]);
              $n=$n+1;

              if ($n==10 && $i<=>($elementos-1)){
                print "\n";
                printf ("VELF%16s");
                $n=0;
               }
           }    
           print "\n"; 
     printf("SPNT%21s\n",$cdp);
     } 
    $elementos=0;
    }
   
   if($linea =~ m/^[0-9]/i && $n_linea>1)
     {
     @line_split= split(/\ +/,$linea);
        $i=0;
    	foreach (@line_split) {
        
        @tvpairs[$i+$elementos]=$_;
        $i=$i+1;
        }
    
     $elementos=$elementos+$#line_split;
     }

}
	  
 $n=0;

          printf ("VELF%17s");
	  for ($i=0;$i<$elementos;$i++) {
            printf ("%5d",@tvpairs[$i]);
            $n=$n+1;
             if ($n==10 && $i<=>($elementos-1)){
               print "\n";
               printf ("VELF%16s");
               $n=0;
               }
           }    
           print "\n";
 
close out;
