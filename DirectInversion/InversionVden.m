%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% How can you invert for the model, given the Green's functions?          %
%                                                                         %
%  \omega^2 u m +Lu = f                                                   %
%                                                                         %
% => \omega^2 u  m = f - Lu                                               %
%                                                                         %
%   m  \sum_{sources} \sum_\omega  \omega^2 u =                           %
%      \sum_{sources} \sum_\omega       f -Lu                             %
%                                                                         %
%                                                                         %
%  m = { \sum_{sources} \sum_\omega f -Lu }  /                            %
%      { \sum_{sources} \sum_\omega  \omega^2 u}                          %
%                                                                         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% define model
n  = [51 101];
h  = [20 20];
z  = [0:n(1)-1]'*h(1);
x  = [0:n(2)-1]*h(2);
[zz,xx] = ndgrid(z,x);

v0 = 2000*ones(n);
dv = zeros(n);
dv = dv - 100*(abs(zz-200)<=10);
%dv = dv + 200*((zz-400).^2+ (xx-1000).^2<=1000);
dv = dv + 200*(abs(zz+.25*xx-1000)<=10);

d0 = 1*ones(n);
dd = zeros(n);
dd = 200*((zz-400).^2+ (xx-1000).^2<=1000)+d0;

s2 = 1./(v0(:)+dv(:)).^2;
b = 1./(dd(:));


% initialized images
num= zeros(n);
den = zeros(n);

% loop over frequencies
sloc = [1:2:100]; %source locations in x 
ns = size(sloc,2);


Ps = getP(n,2,sloc); % source coordinate injection operator
Q  = speye(ns);
F = Ps'*Q;  % source functions
[Dx,Dz] = DifferenceOperators(h,n);


% obtain matrix for 2 different frequencies
f1 = 1.0;
[A1]  = getArho(f1,s2,b,h,n);
f2 = 1.3;
[A2]  = getArho(f2,s2,b,h,n);

N =prod(n);
U1  = A1\F; % green functions for every source function
U2  = A2\F;
sol = zeros(2*N,1);
for is=[1:ns]
    u1 = U1(:,is); % load green's functions for source is
    u2 = U2(:,is);
    fs = F(:,is);% load source function for position is
      
    omega1 = 2*pi*f1;
    M1 = omega1^2*spdiags(u1,0,N,N);
    omega2 = 2*pi*f1;
    M2 = omega2^2*spdiags(u2,0,N,N);
    
    B1 = -Dx'*spdiags(Dx*u1,0,N,N)-Dz'*spdiags(Dz*u1,0,N,N);
    B2 = -Dx'*spdiags(Dx*u2,0,N,N)-Dz'*spdiags(Dz*u2,0,N,N);
    
    lhs = [M1 B1; M2 B2];
    rhs = [fs; fs];
    minv = lhs\rhs;
end

