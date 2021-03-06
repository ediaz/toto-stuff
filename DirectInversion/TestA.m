%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%      test for propagators                                               %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% define model
n  = [50 101];
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
sloc = [3]; %source locations in x 
ns = size(sloc,2);


f = 10.0;
Ps = getP(n,2,sloc);
Q  = speye(ns);
F = Ps'*Q;
[Acden,Scden] = getA(f,s2,h,n,0);
[Avden,Svden] = getArho(f,s2,b,h,n,0);

Ucden = Acden\F;
Uvden = Avden\F;



figure
imagesc(reshape(real(Ucden(:,1)),n));
title('constant den eq');


figure
imagesc(reshape(real(Uvden(:,1)),n));
title('variable den eq');

figure
imagesc(Svden);
colorbar();

figure
imagesc(Scden);
colorbar();