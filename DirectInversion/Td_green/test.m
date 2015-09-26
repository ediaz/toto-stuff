close all;
fileID = fopen('wave.dat','r');
wave = fread(fileID,'float32');
fclose(fileID)

dims = [251,601,476];

wave = reshape(wave,dims);

fqs = [10:15]; % number of frequencies to invert simultaneously
i=0;
for f=fqs
    gf = extractGreen(wave,f,0.004);
    gf = gf(50:end,:);
    imagesc(real(gf));
    pause(0.5)
    %
    n = size(gf);
    %
    u = gf(:); 
    %
    %
    %
    %
    N = prod(n);
    h  = [4 4];
    z  = [0:n(1)-1]'*h(1);
    x  = [0:n(2)-1]*h(2);
    [zz,xx] = ndgrid(z,x);
    %
    %
    L = LaplacianOperator(h,n);
    %
    if(i==0)
        Lhs = spdiags(u.*(2*pi*f)^2,[0],N,N);
        rhs =  -L*u;
    else
        Lhs = [Lhs;spdiags(u.*(2*pi*f)^2,[0],N,N)];
        a =   -L*u;
        rhs = [rhs;a];    
    end
    i=i+1
end

mest = Lhs\rhs;
mest = reshape(abs((real(mest))),n);
figure
imagesc(1./sqrt(mest(2:end-1,2:end-1)));
colorbar()
%
