#! /usr/bin/perl

#list_of_lines.txt es un archivo de variables
open(files,"list_of_lines.txt");

$file="*JOB     BRAZ_01 AAAA
*CALL   DSIN
LABEL   11_post-stack mig post proceso
PKEYLST
BBBB    CCCC
*CALL   GOUT
TAPEOPT /tapefile=\"SEGY/AAAA.sgy\"
*END
";

while ( $linea= <files>)
{

@values = split(/\ +/, $linea); # separa los campos de la linea en el arreglo @values
                                # si se quiere acceder al campo 2 escribir @values[1]
$file=~ s/AAAA/$line_name/g;  #sustitucion de variables AAAA, BBBB, CCCC con los campos
$file=~ s/BBBB/$cdp_in/g;     #del paso anterior
$file=~ s/CCCC/$cdp_fin/g;

$shot_in=@values[2];
print $shot_in;

}

close OUT;
close files;

