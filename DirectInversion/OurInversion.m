%
%
%
%
% How can you invert for the model, given the Green's functions?
%
%
%  \omega^2 u m +Lu = f
%
% => \omega^2 u  m = f - Lu
%
%   m  \sum_{sources} \sum_\omega  \omega^2 u = \sum_{sources} \sum_\omega f -Lu 
%
%  
%  m = { \sum_{sources} \sum_\omega f -Lu }  / { \sum_{sources} \sum_\omega  \omega^2 u} 
%  
%





% define model
n  = [51 101];
h  = [20 20];
z  = [0:n(1)-1]'*h(1);
x  = [0:n(2)-1]*h(2);
[zz,xx] = ndgrid(z,x);

v0 = 2000*ones(n);
dv = zeros(n);
dv = dv - 100*(abs(zz-200)<=10);
dv = dv + 200*((zz-400).^2+ (xx-1000).^2<=1000);
dv = dv + 200*(abs(zz+.25*xx-1000)<=10);

% initialized images
num= zeros(n);
den = zeros(n);

% loop over frequencies
sloc = [1:10:101]; %source locations in x 
ns = size(sloc,2);
for f = [1:10]
    % define operators
    Ps = getP(n,2,sloc);
    [A,L]  = getA(f,1./(v0(:) + dv(:)).^2,h,n); % the getA optionally 
                                                % outputs the Laplacian op
    Q  = speye(ns);

    F = Ps'*Q;  % source functions
    U  = A\(F); % green functions for every source function
    for is=[1:ns]
        u = U(:,is);
        fs = F(:,is); 
        num = num + reshape(conj(u).*(fs-L*u),n);
        den = den + reshape(conj(u).*u.*(2*pi*f)^2,n);
    end
end

mest = real(num.*(1./den)) ; 

vest = 1./sqrt(mest(2:50,2:100));

figure
imagesc(mest); colorbar();
title('inverted model with boundary layer (slowness squared)');


figure
imagesc(vest); colorbar();
title('inverted model without boundary layer');

figure
imagesc(reshape(v0+dv,n)); colorbar()
title('true model');

figure
imagesc(reshape(real(u),n));

