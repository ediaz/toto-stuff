clc;clear;
format long e
N = 64; 
k = 10; 

A = zeros(N,N);
for i=1:N 
    A(i,i) = 2 +1/100;
end
for i=1:N-1
    A(i,i+1) = -1;
end
for i=2:N
    A(i,i-1) = -1;
end
% preconditioner: 
M = zeros(N,N);
for i=1:N 
    M(i,i) = 2;
end

for i=1:N-1
    M(i,i+1) = -1;
end

for i=2:N
    M(i,i-1) = -1;
end
% ============================
x_exact = ones(64,1);
b = A*x_exact;
 
x = zeros(64,1); % initial model
r = b;
z = M\r;
p = z;
for n = 1:k
    if r == zeros(64,1)
        x_exact_obtained = x;
        break;
    end
    s = z'*r;   
    alpha = s/(p'*A*p);
    x = x + alpha*p;
    r = r - alpha*A*p;   
    z = M\r;
    beta = z'*r/s;
    p = z + beta*p;   
    fprintf('n=%2d l2=%16.14e\n',n,norm(x_exact - x))
end
