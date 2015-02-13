#! /usr/bin/perl 
###############################################################################
#
# NAME
#    ttr.pl
#
# SYNOPSIS
#    ttr -t table < stdin > stdout
#
# DESCRIPTION
#    Perl script for performing table-driven character conversions.
#
# AUTHOR
#	phousto@cse.dnd.ca
#
# DATE
#    940327
#
##############################################################################
$THISPROG = 'ttr';
$VERSION = '1.0';
$IDENT = sprintf("%s(%s)", $THISPROG, $VERSION);
$USAGE = sprintf("Usage: %s -t table < stdin > stdout", $THISPROG);
$BUFSIZ = 3200;

# Parse cmd line

while($ARGV[0] =~ /^-/) {
	$_ = shift;
	if(/^-[tT]/) {
		$TBLNAME = $ARGV[0];
	}
}

# Open conversion table file

defined $TBLNAME ||
	die $IDENT, ": ERROR - table file name not provided\n$USAGE\n";
open(TABLE, $TBLNAME) ||
	die "$IDENT: ERROR - Can't open table file '$TBLNAME' ($!)\n";

# Slurp up conversion table

while(<TABLE>) {
	if(/^[^#\/]/) {
		chop;
		push(@tbl, split(/\s*\D+\s*/));
	}
}
close(TABLE);

# Filter STDIN to STDOUT
$strl=3200*8;

#while(read(STDIN, $buf, $BUFSIZ) > 0) {
#
#	for($i = 0; $i < length($buf); $i++) {
#		vec($buf, $i, 8) = $tbl[vec($buf, $i, 8)];
#	
##		if($i=80*8*$j){print substr($buf,($j-1)*80*8
#       }
#	
#$x=join '',$x,$buf;
#}

#for ($i=1; $i< 3200 ;$i=$i+80){
#$tmp= substr($x,$i-1,80);
#print "$tmp\n";
#}
## END OF SCRIPT


read(STDIN, $HEADER, $BUFSIZ) ;
	for($i = 0; $i < length($HEADER); $i++) {
		vec($HEADER, $i, 8) = $tbl[vec($HEADER, $i, 8)];
        }
	
$x=join '',$x,$HEADER;


for ($i=1; $i< 3200 ;$i=$i+80){
$tmp= substr($x,$i-1,80);
print "$tmp\n";
}
# END OF SCRIPT


