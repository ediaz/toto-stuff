clear; clc;
format short e
ns = [10,20,40,80];


u = eps/2;
for i=1:4
  N= ns(i);
  A = randn(N);
  [L1 U1] = gauss(A);
  [L2 U2 P] = lu(A);

  fprintf('--------------------\n')
  fprintf('     N=%d      \n',N) 
  mynorm = norm(L1*U1-A,inf)/norm(A,inf)
  mnorm = norm(L2*U2-P*A,inf)/norm(A,inf)
end


fprintf('--------SMAL a11 --------\n')
for i=1:4
  N= ns(i);
  A = randn(N);
  A(1,1) = u;
  [L1 U1] = gauss(A);
  [L2 U2 P] = lu(A);

  fprintf('--------------------\n')
  fprintf('     N=%d      \n',N) 
  mynorm = norm(L1*U1-A,inf)/norm(A,inf)
  mnorm = norm(L2*U2-P*A,inf)/norm(A,inf)
end

