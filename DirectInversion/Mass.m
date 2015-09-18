function [ M ] = Mass(omega,s2,n,h,freeSurface )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

N=prod(n);

w = ones(n) ;
w(:,1) = 0;    %left
w(:,n(2)) = 0; %rigth
w(1,:) = 0;    %uppper 
w(n(1),:) = 0; %bottom

if freeSurface==1
    w(1,:)=1;
end

% Mass matrix including ABC's


v = (1-w);
v(:,[1 end]) = v(:,[1 end])/h(2);
v([1 end],:) = v([1 end],:)/h(1);

M = omega^2*spdiags(w(:).*s2,0,N,N)+1i*omega*spdiags(v(:).*sqrt(s2),0,N,N);
end

