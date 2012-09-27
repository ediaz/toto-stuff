#!/bin/sh

photo=$1
output=$2


convert -resize 600x750\!  -depth 32 $photo resize_$photo; sleep 0.5s; 
convert resize_$photo -colorspace Gray gray_resize_$photo;   sleep 0.5s;
convert -rotate -90   gray_resize_$photo rot_gray_resize_$photo; sleep 0.5s; 
convert -depth 32  -size 750x600 -endian MSB   -define quantum:format=floating-point \
rot_gray_resize_$photo gray:$output
