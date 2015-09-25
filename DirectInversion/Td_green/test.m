close all;
fileID = fopen('wave.dat','r');
wave = fread(fileID,'float32');
fclose(fileID)

dims = [251,601,476];

wave = reshape(wave,dims);


f=10.0;


gf = extractGreen(wave,f,0.004);
gf = gf(50:end,:);
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
A= spdiags(u.*(2*pi*f)^2,[0],N,N);
rhs = -L*u;
[Dx,Dz] =DifferenceOperators(h,n);
mest = A\rhs;
mest = reshape(abs((real(mest))),n);

imagesc(1./sqrt(mest(2:end-1,2:end-1)));
colorbar()
%
