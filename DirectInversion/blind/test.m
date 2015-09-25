


load('G.mat');
gall = gfull;

f=7.0;

gfT = extractGreen(gfull,f,0.002,1);
gf = gfT';

n = size(gf);

u = gf(:); 




N = prod(n);
h  = [4 4];
z  = [0:n(1)-1]'*h(1);
x  = [0:n(2)-1]*h(2);
[zz,xx] = ndgrid(z,x);


L = LaplacianOperator(h,n);

A= spdiags(u.*(2*pi*f)^2,[0],N,N);
rhs = -L*u;

