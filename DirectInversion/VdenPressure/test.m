close all;
fileID = fopen('/Users/ediaz/Dropbox/toto_y_yo/snap_sp.rsf@','r');
wave = fread(fileID,'float32');
fclose(fileID)

dims = [201,501,476,10];

wave = reshape(wave,dims);

%%

h = [0.004,0.004];
n = [201,501];
N = prod(n);
[Dx,Dz] = DifferenceOperators(h,n);

gf = extractGreen(wave,10,0.004,5);
imagesc(real(gf))
freqs = [11:2:15]; % number of frequencies to invert simultaneously
srcs = [1:10];
nw = size(freqs,2);
ns = size(srcs,2);
i=0;
for f = freqs
    % define operators
    for is=[1:ns]
        omega = 2*pi*f;
        u =extractGreen(wave,f,0.004,is);
        u = u(:);
        B = -Dx'*spdiags(Dx*u,0,N,N)-Dz'*spdiags(Dz*u,0,N,N);
        M = omega^2*spdiags(u,0,N,N);
        if(i==0)
            Lhs = [M B];
            rhs = u*0;
            size(rhs);
            size(Lhs)
        else
            X = [M B];
            Lhs = [Lhs;X];
            rhs = [rhs;u*0];    
        end
        i=i+1;
        fprintf('set %d out of %d\n',i+1,nw*ns);
    end
end


wave=zeros(1); 