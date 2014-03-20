clear; clc;
format  long e

A = pascal(12);
x = ones(12,1);
[L U P]=lu(A);

b = A*x;
c=P*b;
y=forwsub(L,c);
x_1 = backsub(U,y)

relative_error = norm(x-x_1,inf)/norm(x,inf)

k = cond(A,inf);
max_err = k*eps/2

