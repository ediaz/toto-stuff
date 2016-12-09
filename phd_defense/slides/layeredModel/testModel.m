%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%      test for propagators                                               %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear; clc; close all;

% define model
n  = [251 601];
h  = [4 4];
z  = [0:n(1)-1]'*h(1);
x  = [0:n(2)-1]*h(2);
[zz,xx] = ndgrid(z,x);

fileid = fopen('model.bin','r')
v = fread(fileid,[251,601],'float');
fclose(fileid);
s2 = 1./(v(:)).^2;



% initialized images
num= zeros(n);
den = zeros(n);

% loop over frequencies
sloc = [300:300];%:n(2)]; %source locations in x 
ns = size(sloc,2);


ns = prod(size(sloc));
 
N = prod(n);
ns = size(sloc,2);
fqs = [10:10:10];
nw = size(fqs,2);
i=0;
nbox = n;
NB = prod(nbox);
Pbox = getP(n,[1:n(1)],[1:n(2)]);

s0sm = smooth(reshape(s2,n),30);
s0sm = reshape(s0sm,n);

for f = fqs
    % define operators
    Ps = getP(n,2,sloc); % source coordinate injection operator
    [A,L]  = getA(f,s2,h,n,1); % the getA optionally 
                                                % outputs the Laplacian op
    Q  = speye(ns);

    F = Ps'*Q;  % source functions
    U  = A\(F); % green functions for every source function
    for is=[1:ns]
        u = U(:,is); % load green's functions for source is
        ubox = u;
        Lbox = LaplacianOperator(h,nbox);
        
        fs = F(:,is);% load source function for position is
        if(i==0)
            Lhs = spdiags(ubox.*(2*pi*f)^2,[0],NB,NB);
            rhs =  -s0sm(:).*ubox.*(2*pi*f)^2-Lbox*ubox;
            size(rhs)
        else
            Lhs = [Lhs;spdiags(ubox.*(2*pi*f)^2,[0],NB,NB)];
            a =   -s0sm(:).*ubox.*(2*pi*f)^2-Lbox*ubox;
            rhs = [rhs;a];    
        end
        i=i+1;
        fprintf('set %d out of %d\n',i+1,nw*ns);
    end
end
%%


[A,L]  = getA(f,s2(:),h,n,1); % the getA optionally 
                                                % outputs the Laplacian op
Q  = speye(ns);

F = Ps'*Q;  % source functions
U  = A\(F); % green functions for every source function
u = U(:,1); % load green's functions for source is

Fs0sm = fopen('wave.bin','w');
fwrite(Fs0sm,full(real(u)),'float32');
fclose(Fs0sm);

%Invert by least squares:
mest = reshape(real(Lhs\rhs),nbox);
vinv = 1./sqrt(abs(mest));

xrange = [1 (n(2))]*h(2);
zrange = [1 (n(1))]*h(1);
 
figure
% fix boundary:

imagesc(xrange,zrange,mest(3:n(1)-2,3:n(2)-2));
title('inverted model');
daspect([2 n(2)/n(1) 1]);
colorbar()

figure
s2 = reshape(s2,n);
s0sm = reshape(s0sm,n);
imagesc(xrange,zrange,s2(3:n(1)-2,3:n(2)-2)-s0sm(3:n(1)-2,3:n(2)-2));
title('inverted model');
daspect([2 n(2)/n(1) 1]);
colorbar()

Fs2 = fopen('s2.bin','w');
fwrite(Fs2,s2,'float32');
fclose(Fs2);


Fs0sm = fopen('s0sm.bin','w');
fwrite(Fs0sm,s0sm,'float32');
fclose(Fs0sm);


Fs0sm = fopen('ds.bin','w');
fwrite(Fs0sm,full(mest),'float32');
fclose(Fs0sm);

% print('Fig/minv','-depsc');
% 
% figure
% imagesc(xrange,zrange,v(3:n(1)-2,3:n(2)-2));
% title('true model');
% daspect([2 n(2)/n(1) 1]);
% caxis( [1900 2500] )
% 
% colorbar()
% print('Fig/m','-depsc');
% 
% 
% u = real(reshape(u,n));
% 
% figure
% imagesc(xrange,zrange,u(3:n(1)-2,3:n(2)-2)*1000);
% title('wavefield xs=300');
% daspect([2 n(2)/n(1) 1]);
% colorbar()
% 
% print('-depsc','Fig/u');










