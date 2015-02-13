if [ "a$1" == "a" -a "a$2" == "a" ] ; then
    echo "Usage: $0 input.H  surface.pick  a "
    echo "                               "
    echo "a=-1 means remove surface      "
    echo "a=1  means apply surface       "
    echo "	I recommend apply surface with this script because keeps in mind    "
    echo "	the interpolation, window errors                                    "
    echo "	If you apply surface.pick to surface_input.H then                   "

    echo "	datum_surface_input.H== input.H                                     "   
    echo "	if you use VelMerge with surface_input.H and compare with input.H   "
    echo "	you will see a shift between them                                   "
    echo "	This script references the input field from the surface.pick horizon"
    echo "	It was meant to extract the water layer from a velocity file but  can"
    echo "	be used to remove a static for example"
    


   exit 1
fi





velfile=$1  #input
wb=$2
b=$3                        

output=surface_$velfile

if [ "$b" -ne "-1" -a "$b" -ne "1" ] ; then
	 
 	echo "error third input can be 1 or -1"
	exit 1  

   
fi

function getpar() {
    i=$1
    par=$2
    grep "${par}=" $i | tail -n 1 | sed "s/.*${par}=\([^ \t][^ \t]*\).*/\1/"
}

if [ "$b" -eq "-1"  ] ; then
	
	echo "will remove statics"	
	n1_in=`getpar "$velfile" n1`
	o1_in=`getpar "$velfile" o1`
	d1_in=`getpar "$velfile" d1`


	min_depth=`cat $wb|tail -$lines|head -$top|sort -n -k3|head -1|awk '{print $3}'`
	max_depth=`cat $wb|tail -$lines|head -$top|sort -n -r -k3|head -1|awk '{print $3}'`
	 

	npad=`echo "scale=1; $max_depth / $d1_in "|bc -l |awk '{x=$1; y=int(x);if((x-y)>=0.5){x=y+1};x=int(x);print x }'`
 
	
	Pad beg1=$npad end1=$npad <  $velfile |UnGaper dvh=100 intp_type=2 >tmp1_$velfile 
        diff=`echo "scale=1; $npad*$d1_in - $max_depth + $d1_in"|bc -l`

	################################################################################################
	o1_tmp1=`getpar "tmp1_$velfile" o1`
	cat $wb |head -13 >neg_$wb
	lines=` wc $wb|awk '{x=$1-13;print x}'`
	top=`echo "$lines"| awk '{t=$1-1;print t}'`
	echo "cat $wb | tail -$lines| awk '{z=-\$3;printf (\"%10.1f%10.1f%10.1f\n\",\$1,\$2,z)}'|head -$top >>neg_$wb "|sh
	 
	###############################################################################################



	mx_dx2=`echo "scale=0; $max_depth * 2"|bc -l`

	VelMerge hz=neg_$wb vel_file=tmp1_$velfile function=1<tmp1_$velfile >tmp2_$velfile


	o1_tmp2=`getpar "tmp2_$velfile" o1`


	wb_line=`echo "scale=0; -1*$o1_tmp2 + $d1_in "|bc -l`

	
	wb_window=`echo "scale=0; ($wb_line-1*$o1_tmp2)/$d1_in "|bc -l`





	################    Output file a=-1  (remove) ######################################

	Window3d <tmp2_$velfile f1=$wb_window >$output  

	 

	grep n1 $output|tail -1|sed -e "s/$wb_line/0/g" >>$output


	#####################################################################################
	echo "Wb_vel_par  shift_to_keep_in_mind  $diff" >>$output
 

	rm -rf tmp?_$velfile* neg_$wb

else
	echo "will apply statics removed by this script"
  	sfht=`grep "Wb_vel_par  shift_to_keep_in_mind" $velfile|awk '{print $3}'`
        
	VelMerge hz=wb.pick vel_value=1499.0 function=1 <$velfile  >tmp1_$velfile 


	grep n1 tmp1_$velfile|tail -1 |sed -e "s/o1=0.0/o1=-$sfht/g">>tmp1_$velfile

        n1=`grep n1 tmp1_$velfile|tail -1 |sed -e "s/=/ /g"|awk '{print $2}'`
        InterpX o1out=0 n1out=$n1 type=1 <tmp1_$velfile>apply_lin_$velfile
    
        rm -rf tmp1_$velfile*  
fi







