#!/usr/bin/perl

$x1=497161;
$y1=827492;

$x2=497184;
$y2=827481;


$coseno=($x1*$x2 +$y1*$y2)/(sqrt($x1*$x1 +$y1*$y1)*sqrt($x2*$x2 +$y2*$y2));
$theta=($coseno);

print "$theta\n";
print "$theta*180/3.141516\n"; 
