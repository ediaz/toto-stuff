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
###############################################################################

$THISPROG = 'segy_scan';
$VERSION = '1.0';
$IDENT = sprintf("%s(%s)", $THISPROG, $VERSION);
$USAGE = sprintf("Usage: %s < stdin > stdout", $THISPROG);

$check=0;

if  (-t STDIN ) {

 print " Error1: NO INPUT SEGY

 NAME
    segy_map.pl

 DESCRIPTION
    Extract and write to a text file x,y,z  headers which after can be mapped 
    using GMT    
    A combination could be for example Inline, Crossline, elevation

 SYNOPSIS
    segy_map.pl [OPTIONS] < segy to scan > map_xyz.txt

 OPTIONS 

    n=[ntraces]       Number of traces to print, if the user wants to print all traces
                       should put a number considerable bigger that expected number of traces
                       let's say 9999999999999999 (Optional, default= 100)

    s=[start_trace]   Number of FIRST trace to print (Optional, default= 1)


    xheader=[start_byte]    Start byte for x header


    xoffset= [bytes2read]        Number of bytes for x header (normaly either 2 or 4)


    yheader=[start_byte]    Start byte for y header


    yoffset=[bytes2read]        Number of bytes for y header (normaly either 2 or 4)


    zheader=[start_byte]    Start byte for y header


    zoffset=[bytes2read]        Number of bytes for z header (normaly either 2 or 4)




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

while($ARGV[0] =~ /^n/ || $ARGV[0] =~ /^s/ || $ARGV[0] =~ /^xheader/ || $ARGV[0] =~ /^xoffset/ || 
      $ARGV[0] =~ /^yheader/ || $ARGV[0] =~ /^yoffset/ || $ARGV[0] =~ /^zheader/ ||
      $ARGV[0] =~ /^zoffset/ ) {

	$_ = shift;
        if(/^n=/) {
    $aux = $_;
    @array= split(/=/,$aux);
    $traces2print=@array[1];
  }
  if(/^s=/) {
    $aux = $_;
    @array= split(/=/,$aux);
    $first2print=@array[1];
  }
  if(/^xheader=/) {
    $aux = $_;
    @array= split(/=/,$aux);
    $x_header =@array[1];
  }
  if(/^xoffset=/) {
		$aux = $_;
    @array= split(/=/,$aux);
    $x_header_offset =@array[1];
  }     
  if(/^yheader=/) {
		$aux = $_;
    @array= split(/=/,$aux);
    $y_header =@array[1];
  }
  if(/^yoffset=/) {
		$aux = $_;
    @array= split(/=/,$aux);
    $y_header_offset =@array[1];
  }           
  if(/^zheader=/) {
    $aux = $_;
    @array= split(/=/,$aux);
    $z_header =@array[1];
  }
  if(/^zoffset=/) {
    $aux = $_;
    @array= split(/=/,$aux);
    $z_header_offset =@array[1];
  }           
}


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

$nsamples=vec($nsamples,0,16);


#====================> BLOCK 2: CONVERTION AND PRINTING OF EBCDIC AND BINARY HDR<================

#                                              EBCDIC BLOCK:
#ebcdic conv
for($i = 0; $i < length($HEADER); $i++) {
	vec($HEADER, $i, 8) = $tbl[vec($HEADER, $i, 8)];
}



######### TRACE HDRS READ BEGINS ###########

$tr=0;
$pr_tr=0; #Printed traces
$file_size= -s STDIN;
$bytes_read=0;

#print "($file_size - $bytes_read)\n";

while( ($file_size - $bytes_read)>0 ) 
{   
  $bytes2read=$x_header-1;
  read(STDIN, $tmp,$bytes2read);   #byte 1- xheader
  read(STDIN, $x,$x_header_offset);
  
  
  $bytes_read=$bytes_read + $bytes2read + $x_header_offset;
  $bytes2read=$y_header-($x_header + $x_header_offset);
  #print "1: $x_header $bytes2read $bytes_read\n";
  
  read(STDIN, $tmp,$bytes2read);
  read(STDIN, $y,$y_header_offset);
  $bytes_read=$bytes_read + $bytes2read + $y_header_offset;
  $bytes2read=$z_header-($y_header + $y_header_offset);
  
  $a=($y_header + $y_header_offset);
  
  read(STDIN, $tmp,$bytes2read);
  read(STDIN, $z,$z_header_offset);
  $bytes_read=$bytes_read + $bytes2read + $z_header_offset;
  #print "3: $z_header $bytes2read $bytes_read\n";
  
  $bytes2read=240-($z_header + $z_header_offset)+1;
  read(STDIN, $tmp,$bytes2read);
  $bytes_read=$bytes_read + $bytes2read;

  if($x_header_offset == 2){
    $xp=unpack("s>",$x);
  }else{
    $xp=unpack("i>",$x);
  }

  if($y_header_offset == 2){
    $yp=unpack("s>",$y);
  }elsif(y_header_offset == 4){
    $yp=unpack("i>",$y);
  }else{

  }

  if($z_header_offset == 2){
    $zp=unpack("s>",$z);
  }else{
    $zp=unpack("i>",$z);
  }

  $tr=$tr+1;
  $print=0;
  if($tr+1>=$first2print && $pr_tr < $traces2print ){
    printf("%15d%15d%15d\n",$xp,$yp,$zp);
    $pr_tr=$pr_tr+1;
    $print=1;
  }
  
  #SEISMIC SAMPLES READING
  read(STDIN, $tmp,$nsamples*4);
  $bytes_read=$bytes_read +$nsamples*4;		
  last  if($pr_tr==$traces2print );
}

if($pr_tr==0){ 
  print "\n\n\nDid not print trace headers because the optional variable \"start_trace\" is greater than \"TOTAL NUMBER OF TRACES\" of the SEGY FILE!\n\n";
}

