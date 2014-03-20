function [A,b] = problem1(m,n); 

pts = linspace(0,1,m)';

A = zeros(m,n);
b = zeros(m,1);
A(:,1) = 1;
for j=2:n
  A(:,j) = pts.*A(:,j-1);
end

for i=1:m
  b(i) = exp(sin(4*pts(i)));
end
