#! /usr/bin/perl 






$XO=365788.4;
$YO=9666050.8;
$AZ=-143.56;
$Left_handed=1;  #   //1 for left handed -1 for right_handed
$D_inline=30.0;
$D_xline=30.0;
$Inline_orig=2;
$Inline_grid_inc=1;
$Xline_orig=1;
$Xline_grid_inc=1;
$PI=3.141592653589793;

#//For this job, shot and receiver coordinates and headers will be calculated in
#//local coordinates and rotated at the end
#// The suffix LX is local coord X
#// The suffix LY is local coord Y
#// The suffix LOCX is actual coord X
#// The suffix LOCY is actual coord Y
#if ($#ARGV != 0) {
# print "INPUT PARAMETER ERROR!

#       Usage:
#       ./text_coordinates_bin.pl text_file
#
#       Text file is a horizon with columns x,y,z

#       \n";
# exit;
#}


#$entrada=$ARGV[0];

#open(par,"$entrada");

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
		$KCDPY= int(($RCDPYL-0.5*$D_inline)/$D_inline)+ $Inline_orig;
		$KCDPX= int(($RCDPXL-0.5*$D_xline )/$D_xline)+ $Xline_orig;
		        
		$RCDPXLO= ($RCDPX-$XO);
		$RCDPYLO= ($RCDPY-$YO);


		printf("%10s%10s%15.2f\n",$KCDPX,$KCDPY,$Z);
#		printf("%15.2f%15.2f%15.2f%15.2f\n",$RCDPXL,$RCDPYL,$RCDPXLO,$RCDPYLO);
		
	}
	else
	{
		print $_;
	}	



}
