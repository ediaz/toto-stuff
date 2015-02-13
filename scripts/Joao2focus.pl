#!/usr/bin/perl 

#FORMATO VEL FOCUS A FORMATO VEL ANP
#En la linea de comando escribir:
#focus2anp.pl [nombre de archivo de velocidades focus] [nombre de la linea]
#el nombre de la linea no puede tener mas de 11 caracteres
#para salvarlo en un archivo escribir:
#focus2anp.pl [nombre del archivo de velocidades focus] [nombre de la linea] > [nombre de archivo de salida]


$entrada="$ARGV[0]";

open (FILE, $entrada);

$n_linea=0;
$sample=1;
$cdp_ant=-9999999;
$cdp=-99999999;
$h=0;
$elementos=0;
$imprimibles=0;
while ($linea=<FILE>)
{
$n_linea=$n_linea+1;

   if($linea =~ m/^cdp/i) #lineas que comienzan con handvel
    {
    $h=$h+1;
    @cdp_split= split(/\ +/,$linea);
    $cdp_ant=$cdp;
    
    $cdp=@cdp_split[1];
    $i=0;
    $n=0;
   print "\n";
    printf("HANDVEL\t%d\n",$cdp);
      
    $elementos=0;
    $elementos_linea=0;
    }
   
   if($linea =~ m/^cdp/i)
        {
        }
        else
     {
       
      $elementos=$elementos+2;
      @cdp_split= split(/\t/,$linea);
      $vel=substr(@cdp_split[2],0,6); 
      $prof=substr(@cdp_split[1],0,6);    
      
      printf("%d\t%d\t",$prof,$vel);        
      $elementos=$elementos+2;
      $elementos_linea=$elementos_linea+2;
      if($elementos_linea==8){
       $elementos_linea=0;
       print "\n";
      }

}
         

}
	  
 $n=0;
close out;
