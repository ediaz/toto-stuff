clear; clc;
format long e
m = 6;

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
Q_exact = zeros(m);

diary problem1.txt

for j=1:m
  lambda(j) = 4*(sin(j*pi/(2*(m+1))))^2;
  for k=1:m
    Q_exact(k,j) = (sqrt(2/(m+1)))*sin(k*j*pi/(m+1));
  end
end

Q_inv_iteration = zeros(m);

w = zeros(m,1);
sum = zeros(m,1);
diff = zeros(m,1);

v0 = rand(m,1);
v0 = v0/norm(v0);
for j=1:m
  
  w = (T-lambda(j)*eye(m))\v0;
  v = w/norm(w,2);
  
  Q_inv_iteration(:,j) = v;
  diff(j) = norm(Q_exact(:,j)-Q_inv_iteration(:,j));
  sum(j) = norm(Q_exact(:,j)+Q_inv_iteration(:,j));
end

diff = diff
sum = sum

diary off
