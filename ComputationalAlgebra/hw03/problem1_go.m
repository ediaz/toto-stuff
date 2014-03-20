
format long e;
clear;
clc; 

% m,n=15,10
m=100;
n=15;
xn = 2.006787453080206e03;
[A,b]=problem1(m,n);
%x_chol = prob_1bi(A,b)
fprintf('--------------------\n')
fprintf('    m=100,n=15      \n')
x_qr = prob_1bii(A,b);
x_svd = prob_1biii(A,b);
x_mlab = A\b;
fprintf('QR solution x(n)\n')
x_qr(n)
fprintf('SVD solution x(n)\n')
x_svd(n)
fprintf('Matlab solution x(n)\n')
x_mlab(n)

m=10;
n=5;
xn = 2.006787453080206e03;
[A,b]=problem1(m,n);
x_chol = prob_1bi(A,b);
x_qr = prob_1bii(A,b);
x_svd = prob_1biii(A,b);
x_mlab = A\b;

fprintf('--------------------\n')
fprintf('    m=10,n=5      \n')
fprintf('Cholesky solution x(n)\n')
x_chol(n)
fprintf('QR solution x(n)\n')
x_qr(n)




