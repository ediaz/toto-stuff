function [hatQ,hatR] = mgs(A) 
%
% This function does the modified reduced QR
% factorization of A using the Gram-Schmidt
% 
  [m,n] = size(A);
  hatR = zeros(n,n);
  hatQ = zeros(m,n);
  for j = 1:n
    vj = A(:,j);

    for i = 1:j-1
      hatR(i,j) = hatQ(:,i)'*vj;
      vj = vj - hatR(i,j)*hatQ(:,i);  
    end
    hatR(j,j) = norm(vj);
    hatQ(:,j) = vj/hatR(j,j);
  end
end
