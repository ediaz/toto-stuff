#! /usr/bin/perl -w
use Math::Interpolate;

$loc0=0;
$tv_points=0;
$nloc=0;
$loc=0;
while (<STDIN>) 
{
s/RT.VEL/      /g;
@array_split=split(/[ \t]+/,$_);
$array_split[$#array_split]=~ s/
//g;
	if($_  =~ m/^\//i || $#array_split==0){
		 
	}	
	else
	{
		if(/CDP/){
			
			$nloc=$nloc+1;		
			if($nloc==1){$loc0=$loc}

			$loc=$array_split[2];
			$tv_pairs=0;

			if($loc !=$loc0){			
#			print $tv_points;
			print "CDP  $loc \n";
				for($i=1; $i<$tv_points+1 ; $i+=1){						
										
					print "$CDP_VEL[$i][1] $CDP_VEL[$i][2] $i\n";	
				}
			}
		}
		else
		{
			$loc0=$loc;
		 
						
			for ($i=1 ; $i<$#array_split+ 1 ; $i+=2) {
				$tv_pairs=$tv_pairs+1;					
				$CDP_VEL[$tv_pairs][1]=$array_split[$i];
				$CDP_VEL[$tv_pairs][2]=$array_split[$i+1];
				#print "$CDP_VEL[$tv_pairs_line][1] $CDP_VEL[$tv_pairs_line][2] $tv_pairs\n";
			}
			 
			$tv_points=$tv_pairs;
		}	
	}


}
