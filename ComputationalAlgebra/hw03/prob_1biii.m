function [x_svd] = prob_1biii(A,b)
[U,S,V] = svd(A,0);

[m,n] = size(A);
r = size(S,1);
y = zeros(n,1);
c = U'*b;
for i=1:r
  y(i) = c(i)/S(i,i);
end
x_svd = V*y;
