#! /usr/bin/perl 



####MODIFY GRID INFORMATION:  #######
#$XO=365788.4;
#$YO=9666050.8;
#$AZ=-143.56;
#$Left_handed=1;  #   //1 for left handed -1 for right_handed
#$D_inline=30.0;
#$D_xline=30.0;
#$Inline_orig=2;
#$Inline_grid_inc=1;
#$Xline_orig=1;
#$Xline_grid_inc=1;
#######################################


$PI=3.141592653589793;

if  (-t STDIN ) {

 print " Error1: NO INPUT xyz text file

 NAME
    XYZ_text_file_horizon_binning.pl

 DESCRIPTION
    This perl script \"bins\" a three column x y z file following the input
    grid information from the user
    
 SYNOPSIS
    ./XYZ_text_file_horizon_binning.pl  [parameters] < xyz.txt > Inline_Xline_Z.txt

    Example:   
    ./XYZ_text_file_horizon_binning.pl XO=300000.4 \\
    YO=700000.8 \\
    AZ=-170.56 \\
    Left_handed=1 \\
    D_inline=30.0 \\
    D_xline=30.0 \\
    Inline_orig=2 \\
    Inline_grid_inc=1 \\
    Xline_orig=1 \\
    Xline_grid_inc=1 < xyz.txt > Inline_Xline_Z.txt



     You can insert the parameters as input in the command line or write them in a paramter file 


     If the input is X       Y     Z (being Z any attribute), the utput will be 
                     Inline  Xline Z

     The input file can have any text header comment, it won't affect the process


 
 PARAMETERS 
   
     There is no default values, therefore all parameters must be suplied

    XO       []    X coordinate origin
    YO       []    Y coordinate origin
    AZ       []    Azimuth of the grid (this is the angle needed for rotate the inline axis to the positive Y),
                   Mesured counter clock-wise and in degrees

		                      +y
                                        |		
          				|                                        
          				|
                                        |		
          				|
                                        |  In this example the angle would be aprox AZ=~ -140 degrees
          				|
          --------------------------------------------------------- +x
                                       /|\\
          			      / | \\
                                     /  |  \\
          	                    /   |   \\
                                 Xline  |   Inlines
          				|


    Left_handed [] 1 for left handed, -1 for right
    D_inline    []  Distance between inlines
    D_xline     []  Distance between xlines
    Inline_orig []  First inline
    Inline_grid_inc [] Inline increment
    Xline_orig  []  First xline
    Xline_grid_inc  [] Xline increment
    par=parfile [] Parameter file with all previous stated parameters





 AUTHOR
        Esteban Diaz Pantin  

 DATE
    2010.5
 
 Version
    1.0 \n\n";
exit;
}


####MODIFY GRID INFORMATION:  #######
#$XO=365788.4;
#$YO=9666050.8;
#$AZ=-143.56;
#$Left_handed=1;  #   //1 for left handed -1 for right_handed
#$D_inline=30.0;
#$D_xline=30.0;
#$Inline_orig=2;
#$Inline_grid_inc=1;
#$Xline_orig=1;
#$Xline_grid_inc=1;


while($ARGV[0] =~ /^XO/ || $ARGV[0] =~ /^YO/ || $ARGV[0] =~ /^AZ/ || $ARGV[0] =~ /^AZ/ || 
      $ARGV[0] =~ /^Left_handed/ || $ARGV[0] =~ /^D_inline/ || $ARGV[0] =~ /^D_xline/ ||
      $ARGV[0] =~ /^Inline_grid_inc/ ||$ARGV[0] =~ /^Xline_grid_inc/ || $ARGV[0] =~ /^Inline_orig/ ||
      $ARGV[0] =~ /^Xline_orig/ || $ARGV[0] =~ /^par/ ) {

	$_ = shift;
        if(/^XO=/) {

		$aux = $_;
                @array= split(/=/,$aux);
                $XO=@array[1];
                
        }
        if(/^YO=/) {

		$aux = $_;
                @array= split(/=/,$aux);
                $YO=@array[1];
                
        }
        if(/^AZ=/) {

		$aux = $_;
                @array= split(/=/,$aux);
                $AZ=@array[1];
                
        }
        if(/^Left_handed=/) {

		$aux = $_;
                @array= split(/=/,$aux);
                $Left_handed=@array[1];
                
        }
        if(/^D_inline=/) {

		$aux = $_;
                @array= split(/=/,$aux);
                $D_inline=@array[1];
                
        }
        if(/^D_xline=/) {

		$aux = $_;
                @array= split(/=/,$aux);
                $D_xline=@array[1];
                
        }
        if(/^Inline_orig=/) {

		$aux = $_;
                @array= split(/=/,$aux);
                $Inline_orig=@array[1];
                
        }
        if(/^Inline_grid_inc=/) {

		$aux = $_;
                @array= split(/=/,$aux);
                $Inline_grid_inc=@array[1];
                
        }
        if(/^Xline_orig=/) {

		$aux = $_;
                @array= split(/=/,$aux);
                $Xline_orig=@array[1];
                
        }
        if(/^Xline_grid_inc=/) {

		$aux = $_;
                @array= split(/=/,$aux);
                $Xline_grid_inc=@array[1];
                
        }

        if(/^par/) {
  		 
                $aux = $_;
                @array= split(/=/,$aux);
                $par=@array[1];

        }
           
}

if($par){

	#if parameter file was given as input then delete all previously set variables
	undef *Xline_grid_inc;
	undef *Inline_grid_inc;
	undef *Xline_orig;
	undef *Inline_orig;
	undef *D_xline;
	undef *D_inline;
	undef *Left_handed;
	undef *AZ;
	undef *YO;
	undef *XO;
	
	open (FILE, "<$par");

	while (<FILE>){
		if( /XO/ ){
			@array= split(/=/,$_);
	                $XO=@array[1];		
		}
		if( /YO/ ){
			@array= split(/=/,$_);
	                $YO=@array[1];	
		}		
		if( /AZ/ ){
			@array= split(/=/,$_);
	                $AZ=@array[1];	
		}		
		if( /Left_handed/ ){
			@array= split(/=/,$_);
	                $Left_handed=@array[1];	
		}		
		if( /D_inline/ ){
			@array= split(/=/,$_);
	                $D_inline=@array[1];	
		}
		if( /D_xline/ ){
			@array= split(/=/,$_);
	                $D_xline=@array[1];	
		}		
		if( /Inline_orig/ ){
			@array= split(/=/,$_);
	                $Inline_orig=@array[1];	
		}
		if( /Xline_orig/ ){
			@array= split(/=/,$_);
	                $Xline_orig=@array[1];	
		}
		if( /Xline_grid_inc/ ){
			@array= split(/=/,$_);
	                $Xline_grid_inc=@array[1];	
		}
		if( /Inline_grid_inc/ ){
			@array= split(/=/,$_);
	                $Inline_grid_inc=@array[1];	
		}
	
	}
		
}

$XO=~s/
//g;
$YO=~s/
//g;
$AZ=~s/
//g;
$Left_handed=~s/
//g;
$D_inline=~s/
//g;
$D_xline=~s/
//g;
$Inline_orig=~s/
//g;
$Xline_orig=~s/
//g;
$Inline_grid_inc=~s/
//g;
$Xline_grid_inc=~s/
//g;









print "\#XO: $XO\n";
print "\#YO: $YO\n";
print "\#AZ: $AZ\n";
print "\#Left_handed: $Left_handed\n";
print "\#D_inline: $D_inline\n";
print "\#D_xline: $D_xline\n";
print "\#Inline orig: $Inline_orig\n";
print "\#Inline inc $Inline_grid_inc\n";
print "\#Xline orig: $Xline_orig\n";
print "\#Xline inc: $Xline_grid_inc\n";

################   Testing for missing parameters  #################################
if  ($XO eq '' ) { print "error 2 , need XO\n";exit}
if  ($YO eq '' ) { print "error 3 , need YO\n";exit}
if  ($AZ eq '' ) { print "error 4 , need AZ\n";exit}
if  ($Left_handed eq '' ) { print "error 5 , need Left_handed\n";exit}
if  ($D_inline eq '' ) { print "error 6 , need D_inline\n";exit}
if  ($D_xline eq '' ) { print "error 7 , need D_xline\n";exit}
if  ($Inline_orig eq '' ) { print "error 8 , need Inline_orig\n";exit}
if  ($Inline_grid_inc eq '' ) { print "error 9 , need Inline_grid_inc\n";exit}
if  ($Xline_orig eq '' ) { print "error 10 , need Xline_orig\n";exit}
if  ($Xline_grid_inc eq '' ) { print "error 11 , need Xline_grid_inc\n";exit}
#####################################################################################
 




$RCDPX=0;
$RCDPY=0;
while (<STDIN>) 
{

	
@array_split=split(/[ \t]+/,$_);

$test2=@array_split[1];			
			

	for ($i=0 ; $i<$#array_split +1 ; $i+=1) {
		$x=@array_split[$i];
			$condition=0;
		if($x=~ m/^[0-9]/i){
			$condition=1;			
			last ;
		}
		if ($condition==0 && $i==$#array_split +1){print "ERROR:cannot find x coord field";exit}
	}



	if($test2=~ m/^[0-9]/i){

		$RCDPX=@array_split[$i];	
		$RCDPY=@array_split[$i+1];	
		$Z=@array_split[$i+2];	
		$RCDPXL = ($RCDPY-$YO)*sin($AZ*$PI/180) + ($RCDPX-$XO)*cos($AZ*$PI/180);
		$RCDPYL = ($RCDPY-$YO)*cos($AZ*$PI/180) - ($RCDPX-$XO)*sin($AZ*$PI/180);
		$KCDPY= int(($RCDPYL-0.5*$D_inline)/$D_inline)*$Inline_grid_inc+ $Inline_orig;
		$KCDPX= int(($RCDPXL-0.5*$D_xline )/$D_xline)*$Xline_grid_inc+ $Xline_orig;
		        
		$RCDPXLO= ($RCDPX-$XO);
		$RCDPYLO= ($RCDPY-$YO);


		printf("%10s%10s%15.2f\n",$KCDPX,$KCDPY,$Z);
#		printf("%15.2f%15.2f%15.2f%15.2f\n",$RCDPXL,$RCDPYL,$RCDPXLO,$RCDPYLO);
		
	}
	else
	{
		print "#$_";
	}	



}
