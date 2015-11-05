close all;
fileID = fopen('/Users/ediaz/Dropbox/toto_y_yo/snap_sp.rsf@','r');
wave = fread(fileID,'float32');
fclose(fileID)

dims = [201,501,476,10];

wave = reshape(wave,dims);
%%

fqs = [21:1:51]; % number of frequencies to invert simultaneously
i=0;
for sou=1:10
    for f=fqs
        %gfT = extractGreen(gfull(:,:,:,:),f,0.002,sou);
        gprfT = extractGreen(wave,f,0.004,sou);
        
        
        n = size(gprfT);
        gpf = gprfT(:);
        N = prod(n);
        h  = [4 4];
        
        [Dx,Dz] = DifferenceOperators(h,n);
        
        L = LaplacianOperator(h,n);
        
        
        M=spdiags(gpf.*(2*pi*f)^2,0,N,N);
        B=-Dx'*spdiags(Dx*gpf,0,N,N)-Dz'*spdiags(Dz*gpf,0,N,N);
        
  
        
        
        if(i==0)
            Lhs = [M B];
            rhs = 0.*gpf;
            size(rhs);
            size(Lhs);
        else
            
            X = [M B];
            Lhs = [Lhs;X];
            rhs = [rhs;0.*gpf];    
        end
    end
end        
        
%%
  %      [m1 m2]=Lhs\rhs;
        
        
        mest    = real(Lhs\rhs);
ests2b  = reshape(mest(1   :N ),n);
ests2b = ests2b(2:n(1)-1,2:n(2)-1);

estb    = reshape(mest(N+1:2*N),n);
estb = estb(2:n(1)-1,2:n(2)-1);

figure
imagesc(ests2b);
title('estimated s2b');

figure
imagesc(estb);
title('estimated b');

figure
imagesc(1./sqrt(ests2b./estb));
colorbar();
title('estimated velocity');

figure
imagesc(1./estb);
colorbar();
title('estimated density');
