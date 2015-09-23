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

figure 
imagesc(v);

colorbar();

ns = prod(size(sloc));
 
N = prod(n);
ns = size(sloc,2);
fqs = [10:10:10];
nw = size(fqs,2);
i=0;
nbox = [41,161];
NB = prod(nbox);
Pbox = getP(n,[90:130],[180:340]);

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
        ubox = Pbox*u;
        Lbox = LaplacianOperator(h,nbox);
        
        fs = F(:,is);% load source function for position is
        if(i==0)
            Lhs = spdiags(ubox.*(2*pi*f)^2,[0],NB,NB);
            rhs =  -Lbox*ubox;
            size(rhs)
        else
            Lhs = [Lhs;spdiags(ubox.*(2*pi*f)^2,[0],NB,NB)];
            a =   -Lbox*ubox;
            rhs = [rhs;a];    
        end
        i=i+1;
        fprintf('set %d out of %d\n',i+1,nw*ns);
    end
end

%Invert by least squares:
mest = reshape(real(Lhs\rhs),nbox);
figure
imagesc(1./sqrt(mest(2:nbox(1)-1,2:nbox(2)-1)));
colorbar()



























% 
% 
% figure
% imagesc(reshape(real(Ucden(:,1)),n));
% title('constant den eq');
% 
% 
% figure
% imagesc(reshape(real(Uvden(:,1)),n));
% title('variable den eq');
% 
% figure
% imagesc(Svden);
% colorbar();
% 
% figure
% imagesc(Scden);
% colorbar();