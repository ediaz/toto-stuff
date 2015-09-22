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

clear; clc;

% define model
n  = [51 101];
h  = [20 20];
z  = [0:n(1)-1]'*h(1);
x  = [0:n(2)-1]*h(2);
[zz,xx] = ndgrid(z,x);

v0 = 2000*ones(n);
dv = zeros(n);
dv = dv - 100*(abs(zz-200)<=10);
dv = dv + 200*(abs(zz+.25*xx-1000)<=10);

d0 = 1*ones(n);
dd = zeros(n);
dd = 2*((zz-400).^2+ (xx-1000).^2<=1000)+d0;

s2 = 1./(v0(:)+dv(:)).^2;
b = 1./(dd(:));

% loop over frequencies
sloc = [1:2:100]; %source locations in x 
ns = size(sloc,2);

Ps = getP(n,2,sloc); % source coordinate injection operator
Q  = speye(ns);
F = Ps'*Q;  % source functions
[Dx,Dz] = DifferenceOperators(h,n);

N =prod(n);
sol = zeros(2*N,1);
Ps = getP(n,2,sloc); % source coordinate injection operator
Q  = speye(ns);


freqs = [10:10:10];
nw = size(freqs,2);
i=0;
for f = freqs
    % define operators
    A = getArho(f,s2,b,h,n);
    F = Ps'*Q;  % source functions
    U  = A\(F); % green functions for every source function
    for is=[1:ns]
        omega = 2*pi*f;
        u = U(:,is); % load green's functions for source is
        fs = F(:,is);% load source function for position is
        B = -Dx'*spdiags(Dx*u,0,N,N)-Dz'*spdiags(Dz*u,0,N,N);
        M = omega^2*spdiags(u,0,N,N);
        if(i==0)
            Lhs = [M B];
            rhs = fs;
            size(rhs);
            size(Lhs)
        else
            X = [M B];
            Lhs = [Lhs;X];
            rhs = [rhs;fs];    
        end
        i=i+1;
        fprintf('set %d out of %d\n',i+1,nw*ns);
    end n   nn n
end

mest    = real(Lhs\rhs);
ests2b  = reshape(mest(1   :N ),n);
ests2b = ests2b(2:n(1)-1,2:n(2)-1);

estb    = reshape(mest(N+1:2*N),n);
estb = estb(2:n(1)-1,2:n(2)-1);

figure
imagesc(ests2b);
title('estimated s2b');

figure
imagesc(estb);
title('estimated b');

figure
imagesc(1./sqrt(ests2b./estb));
colorbar();
title('estimated velocity');

figure
imagesc(1./estb);
colorbar();
title('estimated density');


figure
imagesc(reshape(s2.*b,n));

figure
imagesc(reshape(b,n));
