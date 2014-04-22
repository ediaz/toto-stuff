clear; clc;
format long e;
m = 4;
T = zeros(m);
for i=1:m
  T(i,i) = 2; % main diagonal
end
for i=1:m-1
  T(i,i+1) = -1; % upper diag
end
for i=2:m
  T(i,i-1) = -1; %lower diag
end

lambda = zeros(m,1);
for j=1:m
  lambda(j) = 4*(sin(j*pi/(2*(m+1))))^2;
end




i=4
T1 = T;
for k = 1:4
  a_m = T1(i,i);
  a_m1 = T1(i-1,i-1);
  b_m1= T1(i-1,i);

  det = sqrt(((a_m-a_m1)/2)^2+b_m1^2);

  if a_m >= a_m1
    mu = (a_m+a_m1)/2+det;
  else
    mu = (a_m+a_m1)/2-det;
  end
  [q,r] = qr(T1(1:i,1:i)-mu*eye(i));
  T1 = r*q+mu*eye(i);
  disp(T1);
  fprintf('T_{4,4}');
  T1(4,4)-2.618033988749895e+00
end
