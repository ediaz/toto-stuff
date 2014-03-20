function [A] = matrixp3A(epsilon)
  A = [1    1    1;...
       epsilon 0 0;...
       0 epsilon 0;...
       0 0 epsilon]; 
end
