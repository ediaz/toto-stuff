close all;
fileID = fopen('wave.dat','r');
wave = fread(fileID,'float32');
fclose(fileID)

dims = [251,601,476];

wave = reshape(wave,dims);

fqs = [5:4:30]; % number of frequencies to invert simultaneously
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
figure
v = squeeze(1./sqrt(mest(2:end-1,2:end-1)));
v = v(:);
for i=1:size(v(:),1)
    if (v(i)>3500)
        v(i)=3500;
    end
end




imagesc(reshape(v,[(n(1)-2) (n(2)-2)]));
colorbar()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

