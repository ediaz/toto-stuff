function [x] = forwsub(L,d)
[n,m] = size(L);

x = zeros(n,1);
x(1) = d(1)/L(1,1);

for i=2:n
  suma=0;
  for j=1:i-1
    suma = suma +L(i,j)*x(j);
  end
  x(i) = (d(i)-suma)/L(i,i);
end

