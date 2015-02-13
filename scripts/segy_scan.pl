#! /usr/bin/perl 
###############################################################################
#
# NAME
#    segy_scan.pl
#
# SYNOPSIS
#    segy_scan.pl -n ntraces -s start_trace < stdin > stdout
#
# DESCRIPTION
#    Segy scanning based on:      SEG-Y rev 1Release 1.02
#
# AUTHOR
#	ediazpantin@gmail.com
#
# DATE
#    2010.1
#
##############################################################################

$THISPROG = 'segy_scan';
$VERSION = '1.0';
$IDENT = sprintf("%s(%s)", $THISPROG, $VERSION);
$USAGE = sprintf("Usage: %s < stdin > stdout", $THISPROG);

$check=0;

#if (-t STDIN && $#ARGV <= 0 ) {
if  (-t STDIN ) {

 print " Error1: NO INPUT SEGY

 NAME
    segy_scan.pl

 DESCRIPTION
    Perl script for Segy* headers scanning, it prints most important
    segy information on the command window: EBCDIC hdr, Binary hdr 
    and some trace headers

 SYNOPSIS
    segy_scan.pl [OPTIONS] < segy to scan > table with header info

 OPTIONS 

    -n [ntraces]       Number of traces to print, if the user wants to print all traces
                       should put a number considerable bigger that expected number of traces
                       let's say 9999999999999999 (Optional, default= 100)

    -s [start_trace]   Number of FIRST trace to print (Optional, default= 1)

    -f [print_seis]    Flag for print seismic samples, 1 = print, 0 = No print(Optional, default=0)



    *Segy scanning based on:      SEG-Y rev 1Release 1.02
    *http://www.seg.org/SEGportalWEBproject/prod/SEG-Publications/Pub-Technical-Standards/Documents/seg_y_rev1.pdf	

 AUTHOR
        Esteban Diaz Pantin ---> ediazpantin\@gmail.com

 DATE
    2010.1
 
 Version
    2.0 break option removed, number of traces is now reached based upon file size\n\n";
exit;
}

#SETTING DEFAULT OPTIONS
$traces2print=100;
$first2print=1;
$seis_flag =0;

$br=0; #Break flag

while($ARGV[0] =~ /^-/) {
	$_ = shift;
	if(/^-[nN]/) {
		$traces2print = $ARGV[0];
	        $_ = shift; 
	}

        if(/^-[sS]/) {
		$first2print = $ARGV[0];
		$_ = shift;
        }
        
        if(/^-[fF]/) {
		$seis_flag =$ARGV[0];
		$_ = shift;
        }
        
 
       
}
if($seis_flag==1||$seis_flag==0){}else{print "ERROR2: -f option ONLY can be \"0\" or \"1\"\n";exit}




#table-driven conversion between EBCDIC and ASCII characters
@tbl=("000","001","002","003","000","009","000","127","000","000","000","011","012","013","014","015",
"016","017","018","000","000","000","008","023","024","025","000","000","028","029","030","031",
"000","000","000","000","000","010","022","027","000","000","000","000","000","005","006","007",
"000","000","021","000","000","000","000","004","000","000","000","000","019","020","000","026",
"032","000","000","000","000","000","000","000","000","000","091","046","060","040","043","033",
"038","000","000","000","000","000","000","000","000","000","093","036","042","041","059","094",
"045","047","000","000","000","000","000","000","000","000","124","044","037","095","062","063",
"000","000","000","000","000","000","000","000","000","096","058","035","064","039","061","034",
"000","097","098","099","100","101","102","103","104","105","000","000","000","000","000","000",
"000","106","107","108","109","110","111","112","113","114","000","000","000","000","000","000",
"000","126","115","116","117","118","119","120","121","122","000","000","000","000","000","000",
"000","000","000","000","000","000","000","000","000","000","000","000","000","000","000","000",
"123","065","066","067","068","069","070","071","072","073","000","000","000","000","000","000",
"125","074","075","076","077","078","079","080","081","082","000","000","000","000","000","000",
"092","000","083","084","085","086","087","088","089","090","000","000","000","000","000","000",
"048","049","050","051","052","053","054","055","056","057","000","000","000","000","000","000");

#===============================> BLOCK 1: READING <============================================


 	#####################################################
	#         READING EBCDIC HEADER                     #
	##################################################### 

print "\n\n#############################         EBCDIC HEADER      ###########################\n";

read(STDIN, $HEADER, 3200) || die "STDIN DOES NOT EXIST\n";;


        ########################################################################
        #                 BINARY HEADER BLOCK READING                          #
        ########################################################################


read(STDIN, $jobid, 4) ; #Job ID, 								bytes 3201-3204
read(STDIN, $line_num, 4) ; #Line number, 							bytes 3205-3208
read(STDIN, $reel_num, 4) ; #Reel number, 							bytes 3209-3212
read(STDIN, $number_traces_x_ensemble, 2) ; #Number of traces per ensemble MANDATORY PSData, 	bytes 3213-3214
read(STDIN, $aux_traces_x_ensemble,2) ; #Sample interval (microsec) MANDATORY PSData, 		bytes 3215-3216
read(STDIN, $sample_int,2) ; #Sample interval (microsec) MANDATORY, 				bytes 3217-3218
read(STDIN, $orig_sample_int,2) ; #Field rec Sample interval (microsec), 			bytes 3219-3220
read(STDIN, $nsamples,2) ; #Number of samples per data trace MANDATORY,				bytes 3221-3222
read(STDIN, $orig_nsamples,2) ; #Orig number of samples X trace          			bytes 3223-3224
read(STDIN, $data_format,2); #Data format, 1 = 4byte IBM, 2 = 4-byte, two's comlement integer
                             #3 = 2-byte, two's complement integer, 4 = ,etc                    bytes 3225-3226
read(STDIN, $Expected_fold,2);#Ensemble fold (e.g. CMP Fold)                                    bytes 3227-3228


#====================> CHECK POINT BINARY HEADER : OK <==========================

read(STDIN, $Trace_sort_code,2); #Trace sort code, 0: unknown, 1: no sort, 2: CDP sort,...,etc  bytes 3229-3230
read(STDIN, $Vert_sum_code,2); #Vertical sum code, from 2 to 32767                              bytes 3231-3232
read(STDIN, $Sweep_f_start,2); #Sweep freq at start                                             bytes 3233-3234
read(STDIN, $Sweep_f_end,2); #Sweep freq at the end                                             bytes 3235-3236
read(STDIN, $Sweep_length,2); #Sweep length                                                     bytes 3237-3238
read(STDIN, $Sweep_type,2); #Sweep type code                                                    bytes 3239-3240
read(STDIN, $Tr_num_Sweep_ch,2); #Trace number of sweep chan                                    bytes 3241-3242
read(STDIN, $Sweep_tr_taper_l_s,2); #Sweep trace taper length at start                          bytes 3243-3244
read(STDIN, $Sweep_tr_taper_l_e,2); #Sweep trace taper length at end                            bytes 3245-3246
read(STDIN, $Taper_type,2); #Taper type: 1=linear, 2=cos^2, 3=other				bytes 3247-3248
read(STDIN, $Corr_data_traces,2); #Correlated data traces               			bytes 3249-3250
read(STDIN, $Binary_gain_recovered,2); #Binary gain recovered 1=yes 2=no		        bytes 3251-3252
read(STDIN, $Amplitude_recovery_method,2); #Amplitude recovery method 1=none,2=SD,3=AGC,4=other	bytes 3253-3254
read(STDIN, $Measurement_system,2); #Measurement system, 1= Meters 2=Feet                       bytes 3255-3256
read(STDIN, $Impulse_signal_pol,2); #Impulse signal polarity                                    bytes 3257-3258
read(STDIN, $Vib_pol_code,2); #Vibratoty polarity code				        	bytes 3259-3260
read(STDIN, $tmp,240); #Unassigned                                                              bytes 3261-3500
read(STDIN, $SEGY_REV,2); #Segy Rev format							bytes 3501-3502
read(STDIN, $FF_zero,2); #Fixed flag zero							bytes 3503-3504
read(STDIN, $tmp,2); #										bytes 3505-3506
read(STDIN, $tmp,94);#UNASSIGNED								bytes 3507-3600

#====================> CHECK POINT BINARY HEADER : OK <==========================




#====================> BLOCK 2: CONVERTION AND PRINTING OF EBCDIC AND BINARY HDR<================

#                                              EBCDIC BLOCK:
#ebcdic conv
for($i = 0; $i < length($HEADER); $i++) {
	vec($HEADER, $i, 8) = $tbl[vec($HEADER, $i, 8)];
}

#PRINTING 3200 EBCDIC HEADER IN ASCII FORMAT (80 characters per line, 40 lines)	
for ($i=1; $i< 3200 ;$i=$i+80){$tmp= substr($HEADER,$i-1,80);print "$tmp\n"}


#                                            BINARY BLOCK:
print "\n\n#######################BINARY BLOCK HDR INFO: ######################################\n\n";
$sample_int=vec($sample_int,0,16);print "sample interval(microsec) :$sample_int\n";
$nsamples=vec($nsamples,0,16);print "Number of samples :$nsamples\n";
$data_format=vec($data_format,0,16);print "Data format :$data_format\n";
$Expected_fold=vec($Expected_fold,0,16);print "Expected FOLD: $Expected_fold\n"; 





#########                                    TRACE HDRS READ BEGINS                                     ###########

$tr=0;
$pr_tr=0; #Printed traces

print "\n\n######################## TRACE HEADER SCAN #########################################\n\n";
print "\tTRACE\t\tRVR_X\t\tRVR_Y\tRVR_el\t\tSRC_X\t\tSRC_Y\t\tSRC_el\t\tCDP_X\t\tCDP_Y\t\tCDP\t\tOFFSET\t    FFID     CHAN     STN\n";

while(read(STDIN, $trace_SEQ,4)>0 ) #trace seq number within line                                 byte   1-4
{   
read(STDIN, $trace_SEQ,4);#Segy trace seq number                                                 byte   5-8
read(STDIN, $FFID,4);     #Field record FFID                                                     byte   9-12
read(STDIN, $Seq_tr_FFID,4);#Seq trace within FFID                                               byte  13-16
read(STDIN, $Stn,4); #Shot Point number SP                                                        byte  17-20
read(STDIN, $CDP,4); #Ensemble number (CDP,CMP,CIG,etc)                                          byte  21-24
read(STDIN, $Seq_tr_CDP,4); #trace number within CDP                                             byte  25-28
read(STDIN, $TR_ID_code,2); #trace ID code 1=Seismic Data, 2=Dead, 3=dummy , etc                 byte  29-30
read(STDIN, $vertical_SUM_TR,2); #Number of vertical summed traces yielding this trace           byte  31-32
read(STDIN, $horizontal_SUM_TR,2);#Number of horizontal summed traces yielding this trace        byte  33-34 
read(STDIN, $Data_use,2); # 1=Production 2=Test                                                  byte  35-36
read(STDIN, $SOFFSET,4); #Signed OFFSET (not absulute)                                           byte  37-40

###########################################sc1   bytes 41-->68 ##############
read(STDIN, $rvr_elev,4); #Reciver elevation                                                     Byte  41-44 *sc1
read(STDIN, $src_elev,4);  #Source elevation							 byte  45-48 *sc1
read(STDIN, $src_depth,4); #Source depth                                                         byte  49-52 *sc1
read(STDIN, $rvr_datum_elev,4); #Datum elev at reciver						 byte  53-56 *sc1
read(STDIN, $src_datum_elev,4); #Datum elev at source						 byte  57-60 *sc1
read(STDIN, $Wdepth_at_src,4); #Water depth at source						 byte  61-64 *sc1
read(STDIN, $Wdepth_at_rvr,4); #Water depth at reciver						 byte  65-68 *sc1
read(STDIN, $sc1,2); # *sc1: scalar who scale the headers with *sc1                              byte  69-70

########################################## sc2   bytes 73-->88 ##############
read(STDIN, $sc2,2); # *sc2: scalar who scale the headers with *sc2                              byte  71-72

read(STDIN, $SRC_X,4); #Source X coordinate                                                      byte  73-76 *sc2
read(STDIN, $SRC_Y,4); #Source Y coordinate                                                      byte  77-80 *sc2
read(STDIN, $RVR_X,4); #ReciverX coordinate                                                      byte  81-84 *sc2
read(STDIN, $RVR_Y,4); #ReciverY coordinate                                                      byte  85-88 *sc2
read(STDIN, $coord_units,2); #Coordinate units, 1= meters or feet, 2 second of arc, etc          byte  89-90

#Refraction statics related hdrs
read(STDIN, $Weath_vel,2 ); #Weathering velocity                                                 byte  91-92
read(STDIN, $sub_Weath_vel,2); #Sub Weathering velocity                                          byte  93-94

#Static related hdrs
read(STDIN, $Uphole_src,2); #Uphole time at source in miliseconds                                byte  95-96 *sc3
read(STDIN, $Uphole_rvr,2); #Uphole time at receiver in miliseconds                              byte  97-98 *sc3
read(STDIN, $src_statics,2); #Source statics applied in miliseconds                              byte  99-100*sc3
read(STDIN, $rvr_statics,2); #Receiver statics applied in miliseconds                            byte 101-102*sc3
read(STDIN, $total_statics,2); #total statics applied in miliseconds                             byte 103-104*sc3

#Recording related hdrs
read(STDIN, $Lag_timeA,2);# initiation pules that may be recorded on aux trace                   byte 105-106*sc3
read(STDIN, $Lag_timeB,2);#Time in ms between time break (Lag_timeA) and init of source          byte 107-108*sc3
read(STDIN, $Delay_time,2);# Delay rceording time                                                byte 109-110*sc3

#Mute hdrs
read(STDIN, $Start_mute_time,2); #Start mute time                                                byte 111-112*sc3
read(STDIN, $End_mute_time,2); #End mute time                                                    byte 113-114*sc3

#Trace legth hdrs
read(STDIN, $Nsamples_tr,2);#Nsamples for this trace                                             byte 115-116
read(STDIN, $SI_tr,2); #Sample interval for this trace (microseconds)                            byte 117-118

#gain related hdrs
read(STDIN, $Gain_type_of_instrument,2);#gain type of field instrument 1=fixed 2=binary,3=fp,etc byte 119-120
read(STDIN, $instrument_gain_constant,2); #instrument_gain_constant (dB)                         byte 121-122
read(STDIN, $instrument_initial_gain,2); #instrument_initial gain   (dB)                         byte 123-124

#Vibroseis related hdrs
read(STDIN, $correlated,2); #Correlated? 1=yes 2=no                                              byte 125-126
read(STDIN, $Sweep_f_start_tr,2); #Sweep freq start at trace (is also on binary header)          byte 127-128
read(STDIN, $Sweep_f_end_tr,2);   #Sweep freq end at trace (is also on binary header)            byte 129-130
read(STDIN, $Sweep_length_tr,2);#Sweep length at trace   (is also on binary header)              byte 131-132
read(STDIN, $Sweep_type_tr,2); #Sweep type at trace(is also on binary header) 1=lin 2=par 3=exp  byte 133-134
read(STDIN, $Sweep_tr_taper_l_s_tr,2); #Sweep trace taper length at start at trace               byte 135-136
read(STDIN, $Sweep_tr_taper_l_e_tr,2); #Sweep trace taper length at end at trace                 byte 137-138
read(STDIN, $Seep_taper_type_tc,2); #Sweep taper type at trace 1=lin 2=cos^2 3=other             byte 139-140

#Filter related hdrs
read(STDIN, $alias_filter_fr,2); #Alias filter freq (if used)                                    byte 141-142
read(STDIN, $alias_filter_slope,2); #Alias slope (dB/octave)                                     byte 143-144
read(STDIN, $Notch_filter_fr,2); #Notch filter freq (if used)                                    byte 145-146
read(STDIN, $Notch_filter_slope,2); #Notch slope (dB/octave)                                     byte 147-148
read(STDIN, $Low_cut_fr,2); #Low cut freq (if used)                                              byte 149-150
read(STDIN, $High_cut_fr,2); #High cut freq (if used)                                            byte 151-152
read(STDIN, $Low_cut_slope,2);  #Low cut filter slope (dB/octave)                                byte 153-154
read(STDIN, $High_cut_slope,2); #High cut filter slope (dB/octave)                               byte 155-156

#Recording date related hdrs
read(STDIN, $Year,2); #Year in four digits                                                       byte 157-158
read(STDIN, $Julian_day,2); #Julian day                                                          byte 159-160
read(STDIN, $Hour,2); #24 hours based HOUR                                                       byte 161-162
read(STDIN, $Minute,2); #minute of hour                                                          byte 163-164
read(STDIN, $Sec,2); #Second of minute                                                           byte 165-166
read(STDIN, $Time_code,2); #Time code, 1=Local 2=GMT (Greenwich) 3=other 4=UTC(coord.univer.time)byte 167-168
read(STDIN, $Tr_weight,2); #Trace weigth                                                         byte 169-170


#useless geophone related hdrs
read(STDIN, $tmp,10); #USELESS HDRS                                                               byte 171-180

#CDP related hdrs
read(STDIN, $CDP_X,4); #CDP X coordinate                                                          byte 181-184
read(STDIN, $CDP_Y,4); #CDP Y coordinate                                                          byte 185-188
read(STDIN, $Iline,4); #Inline number (3D post stack data)                                        byte 189-192
read(STDIN, $Xline,4); #Crossline number(3D post stack data)                                      byte 193-196
read(STDIN, $SP,4);  #Shot point (Nearest SP to the CDP location)                                 byte 194-200
read(STDIN, $sc_SP,2); #Scalar to be applied to SP, negative divides, positive multiplies         byte 201-202

#Trace measurement related hdrs
read(STDIN, $Tr_measurement,2); #Trace value measurement, 1=pascal,2Volts,3mV,4Amperes,etc        byte 203-204
read(STDIN, $Transduction_Constant_A,4); #Constant to convert between units of measurement A      byte 205-208
read(STDIN, $Transduction_Constant_B,2); #Constant to convert between units of measurement B      byte 209-210
                                         #Trace transduction constant= A*10^B
read(STDIN, $Tr_measurement2,2); #Trase measurement after Transduction Constant applied           byte 211-212


read(STDIN, $Device_id,2); #Device identifier                                                     byte 213-214
read(STDIN, $sc3,2); #Scalar to be applied to headers 95-114 TIME SCALAR                          byte 215-216


#Source measurement related hdrs
read(STDIN, $Src_type_n_orientation,2); #1Vib-verticl,2Vib-Crossline,3VibInline,4ImpulsiveVert... byte 217-218
read(STDIN, $Src_Dir,6); #Energy direction in DEGREES ?                                           byte 219-224
read(STDIN, $Src_Measurement_A,4); #Source effort in physics units,                               byte 225-228
read(STDIN, $Src_Measurement_B,2); #Source effort in physics units: A*10^B                        byte 229-230
read(STDIN, $Src_Measurement_unit,2); #1:Joule, 2:kW,3Pa,4Bar,4Bar-meter(Bar-m),5Newton,6Kg       byte 231-232
read(STDIN, $tmp,8); #Unassigned bytes, if used MUST be specified in EBCDIC hdr                   byte 233-240







###################################################################################################
#BINARY ---- 2 ---> Human readable transformation of table variables:

#Scalar Constants: negatives divides, positives multiplies
$sc1=vec($sc1,0,16);#if($sc1<0){$sc1=-1/$sc1};if($sc1==0){$sc1=1}; #transform for always perform multiplication
$sc2=vec($sc2,0,16);#if($sc2<0){$sc2=-1/$sc2};if($sc2==0){$sc2=1};
$sc3=vec($sc3,0,16);#if($sc3<0){$sc3=-1/$sc3};if($sc3==0){$sc3=1};

#Coordinates transformation
$SRC_X=vec($SRC_X,0,32);#/$sc2;
$SRC_Y=vec($SRC_Y,0,32);#/$sc2;

$RVR_X=vec($RVR_X,0,32);#/$sc2;
$RVR_Y=vec($RVR_Y,0,32);#/$sc2;

$CDP_X=vec($CDP_X,0,32);#/$sc2;
$CDP_Y=vec($CDP_Y,0,32);#/$sc2;



#Elevation transformation
$rvr_elev=vec($rvr_elev,0,32)/$sc1;
$src_elev=vec($src_elev,0,32)/$sc1;

#OFFSET
$SOFFSET=vec($SOFFSET,0,32);
$CDP=vec($CDP,0,32);
$Stn=vec($Stn,0,32);
$FFID=vec($FFID,0,32);
$Seq_tr_FFID=vec($Seq_tr_FFID,0,32);

$tr=$tr+1;

#TRACE\tRVR_X\tRVR_Y\tRVR_el\tSRC_X\tSRC_Y\tSRC_el\tCDP_X\tCDP_Y\tOFFSET\tFFID\tStation\n
$imprime=0;
   if($tr>=$first2print && $pr_tr < $traces2print ){
	printf("tr:%8d\t%8.1f\t%8.1f\t%8.1f\t%8.1f\t%8.1f\t%8.1f\t%8.1f\t%8.1f\t%8.1f\t%8d%8d%8d%8d\n",
	$tr,$RVR_X,$RVR_Y,$rvr_elev,$SRC_X,$SRC_Y,$src_elev,$CDP_X,$CDP_Y,$CDP,
	$SOFFSET,$FFID,$Seq_tr_FFID,$CDP);
        $pr_tr=$pr_tr+1;
	$imprime=1;
    }


#SEISMIC SAMPLES READING
$Nsamples_tr=vec($Nsamples_tr,0,16);
	
	if($seis_flag==1 && $imprime ==1){
	
		for ($i=1; $i<$Nsamples_tr+1; $i=$i+1){
			read(STDIN, $seis,4);
			$seis = unpack("B*", $seis);
			#printf("\t%40.1f  sample:%5d\n",$seis,$i);
			$seis = conv_bit_string_2_ibm32float($seis);		
			printf("\t%40.1f  sample:%5d\n",$seis,$i);
		}
	}
	else {
	
		read(STDIN, $tmp,$Nsamples_tr*4);
		
	}

        last  if($pr_tr==$traces2print )

}

if($pr_tr==0){ print "\n\n\nDid not print trace headers because the optional variable \"start_trace\" is greater than \"TOTAL NUMBER OF TRACES\" of the SEGY FILE!\n\n";
}



print "####################################################################################\n\n";


$file_size= -s STDIN;
print "segy size : $file_size bytes\n";
$traces_from_size=($file_size -3600)/(240+$Nsamples_tr*4);

 

print "TOTAL NUMBER OF TRACES=  $traces_from_size \n";
#USED SUB ROUTINES GRABBED FROM: http://bytes.com/topic/perl/answers/50339-ibm-32-bit-floating-point-conversion


#================================================= =====================
#copied from "Perl Cookbook, recipe 2.4
sub bin2dec {
return unpack("N", pack("B32", substr("0" x 32 . shift, -32)));
}
#================================================= ======================


#================================================= ======================
sub conv_bit_string_2_ibm32float {
$bit_string = shift;

#================================================= ============
#DEBUG (make $bit_string = "11000010101100011001111110101111")
# (1st bit, sign bit = 1 (denotes negative number))
# (next seven bits, exponent = 66)
# (last 24 bits represent fraction = 0.69384283 in decimal)
# (whole bit string as IBM 32 bit FP = -177.62376 decimal)
#================================================= ============

# A value of 1 for sign bit denotes a negative number
# while a value of zero denotes a positive number
#================================================= ===
$first_digit = substr($bit_string, 0, 1);
if ( $first_digit eq "0" ) {
$sign_bit = 1;
} elsif ( $first_digit eq "1" ) {
$sign_bit = -1;
}

$bin_exponent = substr($bit_string, 1, 7);
$exponent = bin2dec($bin_exponent);

#================================================= ====================
#Computing fraction
#The following will do this:
# take last 24 bits of $bit_string such as "101100011001111110101111"
# for each bit starting from left to right, convert to decimal
# i.e. (1 * 2**-1) + (0 * 2**-2) + (1 * 2**-3) + (1 * 2**-4) ...
#================================================= ====================
$bin_fraction = substr($bit_string, 8, 24);
@bit_chars = unpack("A1" x length($bin_fraction), $bin_fraction);

$place_holder = -1;
$fraction = 0;
foreach $bit ( @bit_chars ) {
$fraction += $bit * (2 ** $place_holder);
$place_holder += -1;
}

$ibm_float = ($sign_bit ** 1) * (16 ** ($exponent - 64)) *
($fraction);
return sprintf("%.10f", $ibm_float);
}

__END__


























