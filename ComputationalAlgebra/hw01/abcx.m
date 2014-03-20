function [A,B,C,X] = abc(n)
%
% This function returns 5 matrices A,B,C,X,T
%
%
%

  X = zeros(n,n);
  for i=1:n
    for j=1:n
      X(i,j) = sqrt(2./(n+1))*sin(i*j*pi/(n+1));
    end
  end

  T = zeros(n,n);
  for i=1:n
    T(i,i) = 2.0; %main diagonal
  end
  for i=1:n-1
    T(i+1,i  ) = -1; % lower diagonal
    T(i  ,i+1) = -1; % upper diagonal
  end

  A = X -X';
  B = X*X;
  C = X*T*X;
end 
