clear; clc;
format long e
A = zeros(64,64);
[n,m] = size(A);

for i=1:n
    A(i,i)=i;
end
for i=1:n-1
    A(i,i+1)=1;
end
for i=2:n
    A(i,i-1)=1;
end
eigA = eig(A);
b = ones(n,1);
k=40;
Q = zeros(n, k+1);
Q(:,1) = b/norm(b);
v = A*Q(:,1);
H(1,1) = Q(:,1)'*v;
v = v - H(1,1)*Q(:,1);
H(2,1) = norm(v);
if H(2,1) == 0
    break;
end
H(1,2) = H(2,1);
Q(:,2) = v/(H(2,1));
for j = 2:k
    v = A*Q(:,j);
    H(j-1,j) = Q(:,j-1)'*v;
    v = v - H(j-1,j)*Q(:,j-1);
    H(j,j) = Q(:,j)'*v;
    v = v - H(j,j)*Q(:,j);
    H(j+1,j) = norm(v);  
    if H(j+1,j) == 0
        break;
    end
    Q(:,j+1) = v/(H(j+1,j));  
end
 
T = H(1:k,1:k);
eigT = eig(T);

min_lambda_A = min(eigA)
min_lambda_T = min(eigT)

max_lambda_T = max(eigT) 
max_lambda_A = max(eigA)

min_diff = abs(min_lambda_A - min_lambda_T)
max_diff = abs(max_lambda_A - max_lambda_T)