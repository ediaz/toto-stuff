function [L,U] = gauss(A)

[n,m] = size(A);
U=A;
L=eye(m);

for k=1:m-1
  for i=k+1:m
    multiplier = U(i,k)/U(k,k);
    L(i,k) = multiplier;
    for j=k:m
      U(i,j)=U(i,j)-multiplier*U(k,j);
    end
  end
end
