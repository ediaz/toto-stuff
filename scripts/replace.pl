#!/usr/bin/perl -w

#sintaxis: replace.pl <archivo>
#cambiar el campo AAAA por el BBBB
open (FILE, "$ARGV[0]");

while (<FILE>)
{
#s/AAA/BBB/g
s/292/008/g;
chomp;
push @file, $_;
}

close (FILE);
open (OUT, ">$ARGV[0]");
foreach $x (@file)
{
print OUT "$x\n";
}
