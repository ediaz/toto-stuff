function [x] = backsub(U,d)
[n,m] = size(U);

x = zeros(n,1);
x(n)= d(n)/U(n,n);

for i=n-1:-1:1
  suma = 0;
  for j=i+1:n
    suma = suma+U(i,j)*x(j);
  end
  x(i) = (d(i)-suma)/U(i,i);
end

