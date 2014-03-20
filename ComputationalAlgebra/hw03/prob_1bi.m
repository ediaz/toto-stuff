function [x_chol] = problem_1bi(A,b);
R = chol(A'*A);
c = A'*b;
y = forwsub(R',c);
x_chol = backsub(R,y);
