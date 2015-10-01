close all;

load('gfull.mat');

fqs = [21:4:61]; % number of frequencies to invert simultaneously
i=0;
for sou=1:6
    for f=fqs
        gfT = extractGreen(gfull(:,:,:,:),f,0.002,sou);
        %gfT = One_suo_Green(Green,f,0.004);
        gf=gfT';
        imagesc(real(gf));
        pause(0.5)

        n = size(gf);
        u = gf(:); 
        N = prod(n);
        h  = [4 4];
        L = LaplacianOperator(h,n);
        if(i==0)
            A = spdiags(u.*(2*pi*f)^2,[0],N,N);
            b =  -L*u;
        else
            A = [A;spdiags(u.*(2*pi*f)^2,[0],N,N)];
            a =   -L*u;
            b = [b;a];    
        end
        i=i+1
    end
end


mest = A\b;

%%
[Dx,Dz] = DifferenceOperators(h,n);

xrange = [2 (n(2)-2)]*h(2);
zrange = [2 (n(1)-2)]*h(1);

eps =0.1;
eps1=1E7;

R = eps*Dz;
G2= eps1*Dx;
C = [A;R;G2];
rhs = [b;zeros(N,1);zeros(N,1)];


mestr = C\rhs;
mest  = reshape(abs((real(mest))),n);
mestr = reshape(abs((real(mestr))),n);


figure
%imagesc(mest);
imagesc(xrange,zrange,1./sqrt(mest(2:end-1,2:end-1)));
title('inverted model');
daspect([1 1 1]);
caxis( [1900 2500] )
colorbar()
print('Fig/minv_6','-depsc');


figure;
%imagesc(mestr);
imagesc(xrange,zrange,1./sqrt(mestr(2:end-1,2:end-1)));
title('regularized inverted model');
daspect([1 1 1]);
caxis( [1900 2500] )
colorbar()
print('Fig/minv_reg_6','-depsc');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

