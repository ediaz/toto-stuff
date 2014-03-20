clear;
clc;
format long e;
for i=1:5
  x = (100)^i
  fx1=sqrt(x+1)-sqrt(x)
  fx2=1/(sqrt(x+1)+sqrt(x))
  err = abs(fx2-fx1)/fx2*100
end
