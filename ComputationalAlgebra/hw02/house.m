function [hatQ,hatR] = house(A)
  
  [m,n] = size(A);
  V = zeros(m,n);

  for k=1:n
    x = A(k:m,k);
    v = x + (sign(x(1,1))*norm(x)*eye(m-k+1,1));
    V(k:m,k) = v/norm(v);
    A(k:m,k) = -sign(x(1,1))*norm(x)*eye(m-k+1,1);
    
    for j= k+1:n
      A(k:m,j) = A(k:m,j)-(V(k:m,k)'*A(k:m,j)*2)*V(k:m,k);
    end
  end
  hatR = A(1:n,1:n);
  hatQ = eye(m,n);

  for j=1:n
    for k=j:-1:1
      hatQ(k:m,j) = hatQ(k:m,j)-V(k:m,k)*(V(k:m,k)'*hatQ(k:m,j)*2);  
    end
  end
end
