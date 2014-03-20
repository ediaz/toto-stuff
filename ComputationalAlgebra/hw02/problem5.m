clear all;
clc;
format short e;
diary problem5.txt
m = 30;
n = 20;
pts = linspace(1/n,1,n);
A = zeros(m,n);

A(1,:) =ones(1,20);
for i=2:m
  A(i,:) = pts.*A(i-1,:);
end

[cQ,cR]=mgs(A);
[hQ,hR]=house(A);
[mQ,mR]=qr(A,0);
I = eye(size(mQ'*mQ));

fprintf('   ||A - QR||:\n   mgs = %e   house = %e  matlab = %e \n',...
            norm(A-cQ*cR),norm(A-hQ*hR),norm(A-mQ*mR));
fprintf('   ||I - Q^T*Q||:\n   mgs = %e   house = %e  matlab = %e \n',...
            norm(I-cQ'*cQ),norm(I-hQ'*hQ),norm(I-mQ'*mQ));

diary off;
