function [A,S,M] = getArho(f,s2,b,h,n,freeSurface)
% 5-point discretization of the 2D Helmholtz operator with absorbing
% boundary conditions.
%
% use:
%   A = getA(f,m,h,n);
%
% input:
%   f - frequency [Hz]
%  s2 - slowness squared [m^2/s^2]
%   b - 1/density [kg/m^3]
%   h - [dz,dx] gridspacing in z and x direction [m]
%   n - [nz,nx] number of gridpoints in z and x direction
%
% ouput:
%   A - sparse matrix
%

switch nargin
    case 5
        freeSurface = 1;
end

m1  = s2(:).*b(:);
m2 = b(:);

% angular frequency
omega = 2*pi*f;

% number of gridpoints
N     = prod(n);

% Stiffness matrices
B = spdiags(b,[0],N,N);
[Dx,Dz] = DifferenceOperators(h,n);
S =  -Dz'*B*Dz -Dx'*B*Dx;

% absorving 1 grid point layers (w =1 inside domain, 0 at the boundary layer) 
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
M = omega^2*spdiags(w(:).*m1,0,N,N) + ...
    1i*omega*spdiags(v(:).*sqrt(s2),0,N,N); % boundary condition

% 
A = M + S;




