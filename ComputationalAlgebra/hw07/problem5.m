clc;
clear;
format long e
 
N = 64; 
k = 30; 

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
x_exact = ones(N,1);
b = A*x_exact;
x = zeros(N,1); % initial model
r = b;
p = r;
for n = 1:k
    if r == zeros(64,1)
        x_exact_obtained = x;
        break;
    end    
    z = r'*r;    
    alpha = z/(p'*A*p);
    x = x + alpha*p;    
    r = r - alpha*A*p;  
    beta = r'*r/z;
    p = r + beta*p;    
    if rem(n,5) == 0
        fprintf('n=%2d l2=%16.14e\n',n,norm(x_exact - x))
    end
end
