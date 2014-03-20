function [x_qr] = prob_1bii(A,b)
[Q,R] = qr(A,0); 
c = Q'*b;
x_qr = backsub(R,c);

