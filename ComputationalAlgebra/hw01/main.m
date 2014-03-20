clear;
clc; 

format short e; 

diary output.txt 
for n=2:4 
  fprintf('number of columns n=%d\n',n);
  [A,B,C] = abc(n)
end
diary off
