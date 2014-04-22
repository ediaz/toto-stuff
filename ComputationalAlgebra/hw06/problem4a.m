clear; clc;
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


diary problem4a.txt
T_0 = T;
for k = 1:40
  [Q,R] = qr(T_0);
  T_0 = R*Q;
  if rem(k,10) == 0
    fprintf('i=%d\n',k);
    t = diag(T_0);
    disp(T_0);
    
    fprintf('norm (convergence):');
    disp(norm(sort(t)-lambda)); 
    fprintf('norm (T):');
    disp(norm(T_0-T_0')); 

  end
end
