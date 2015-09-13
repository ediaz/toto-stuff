function [f, g ] = jof(mk,laplac,u,s,omega)
%
% this function computes the gradient and value of
%
% J = || Lm -d || 
%
% where L = diag(omega^2 U)
% and d = f-laplacU
%
% input:
%   mk: current model (column vector)
%   laplac: matrix with space derivatives from the helmotz solver
%   u: green's function for source s 
%   s: source function (same dimension of u)
%   omega: angular frequency
%
% output:
%
%   f: value of objective function for mk model
%   g: gradient of J at model mk
%

n= prod(size(mk));
L = spdiags(omega^2*u,n,n); 
d = s - laplac*u; 

g = L'*(L*m -d); 
f = norm(Lm-d,2); 


end

