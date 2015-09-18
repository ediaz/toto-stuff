function [Dx,Dz] = DifferenceOperators(h,n)



D1 = spdiags(ones(n(1),1)*[1 -1]/h(1),[0 1],n(1),n(1)); 
D1(1,1)=0;

D2 = spdiags(ones(n(2),1)*[1 -1]/h(2),[0 1],n(2),n(2)); 
D2(1,1)=0;

Dz = kron(speye(n(2)),D1');
Dx = kron(D2',speye(n(1)));

end

