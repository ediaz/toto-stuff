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
dv = dv + 200*((zz-400).^2+ (xx-1000).^2<=1000);
dv = dv + 200*(abs(zz+.25*xx-1000)<=10);

% initialized images
num= zeros(n);
den = zeros(n);

% loop over frequencies
sloc = [1:5:101]; %source locations in x 
ns = size(sloc,2);
for f = [1:10]
    % define operators
    Ps = getP(n,2,sloc); % source coordinate injection operator
    [A,L]  = getA(f,1./(v0(:) + dv(:)).^2,h,n,0); % the getA optionally 
                                                % outputs the Laplacian op
    Q  = speye(ns);

    F = Ps'*Q;  % source functions
    U  = A\(F); % green functions for every source function
    u = sum(U,2); % create a normal incidence plane wave
    num = num + reshape(conj(u).*(-L*u),n);
    den = den + reshape(conj(u).*u.*(2*pi*f)^2,n);
    
end


%
% the inverted model is position independent, so actually we don't
%



mest = real(num.*(1./den)) ; 
mest2 = mest;
mest2(2,sloc) = 1/2000^2;

vest = real(1./sqrt(mest(2:50,2:100)));
vest2 = real(1./sqrt(mest2(2:50,2:100)));

figure
imagesc(mest); colorbar();
title('inverted model with boundary layer (slowness squared)');


figure
imagesc(vest); colorbar();
title('inverted model without boundary layer no source');

figure
imagesc(vest2); colorbar();
title('inverted model without boundary layer source filled');

figure
imagesc(reshape(v0+dv,n)); colorbar()
title('true model');

figure
imagesc(reshape(real(u),n));
title('R(u(101,x,f=10)');



figure 
subplot(2,1,1);
imagesc(real(num(2:50,2:100)));
title('Q \sum_\omega -\nabla^2 u');
subplot(2,1,2);
imagesc(den(2:50,2:100));
title('Q \sum_\omega \omega^2 u');



