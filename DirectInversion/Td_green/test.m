close all;
fileID = fopen('wave.dat','r');
wave = fread(fileID,'float32');
fclose(fileID)

dims = [251,601,476];

wave = reshape(wave,dims);

fqs = [11:2:31]; % number of frequencies to invert simultaneously
i=0;
for f=fqs
    gf = extractGreen(wave,f,0.004);
    gf = awgn(gf(50:end,:),-10);
    imagesc(real(gf));
    pause(0.5)
   
    n = size(gf);
    u = gf(:); 
    N = prod(n);
    h  = [4 4];
    z  = [0:n(1)-1]'*h(1);
    x  = [0:n(2)-1]*h(2);
    [zz,xx] = ndgrid(z,x);
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

mest = A\b;
mest = reshape(abs((real(mest))),n);

vinv = 1./sqrt(mest);

%%

xrange = [3 (n(2)-3)]*h(2);
zrange = [3 (n(1)-3)]*h(1);
 
figure
imagesc(xrange,zrange,vinv(4:n(1)-3,4:n(2)-3));
title('inverted model');
daspect([4 n(2)/n(1) 1]);
caxis( [1900 2500] )
colorbar()
print('Fig/minv_10f','-depsc');











