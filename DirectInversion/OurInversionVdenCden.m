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
N = prod(n);
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
dd = 2*((zz-400).^2+ (xx-1000).^2<=1000)+d0;

s2 = 1./(v0(:)+dv(:)).^2;
b = 1./(dd(:));


% initialized images
num= zeros(n);
den = zeros(n);

% loop over frequencies
sloc = [1:5:101]; %source locations in x 
ns = size(sloc,2);
fqs = [10:10];
nw = size(fqs,2);
i=0;
for f = fqs
    % define operators
    Ps = getP(n,2,sloc); % source coordinate injection operator
    Arho =getArho(f,s2,b,h,n,0);
    [A,L]  = getA(f,1./(v0(:) + dv(:)).^2,h,n,0); % the getA optionally 
                                                % outputs the Laplacian op
    Q  = speye(ns);

    F = Ps'*Q;  % source functions
    U  = Arho\(F); % green functions for every source function
    for is=[1:ns]
        u = U(:,is); % load green's functions for source is
        fs = F(:,is);% load source function for position is
        if(i==0)
            Lhs = spdiags(u.*(2*pi*f)^2,[0],N,N);
            rhs =  -L*u;
            size(rhs)
        else
            Lhs = [Lhs;spdiags(u.*(2*pi*f)^2,[0],N,N)];
            a =   -L*u;
            rhs = [rhs;a];    
        end
        i=i+1;
        fprintf('set %d out of %d\n',i+1,nw*ns);
    end
end

%Invert by least squares:
mest = reshape(real(Lhs\rhs),n); 

%
% the inverted model is position independent, so actually we don't
% need to know much about the source (given that is injected in one point).
%



mest2 = mest;
mest2(2,sloc) = 1/2000^2;

vest = 1./sqrt(mest(2:50,2:100));
vest2 = 1./sqrt(mest2(2:50,2:100));

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