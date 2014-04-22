clc;clear;

T = [2 -1 0 0 0 0; -1 2 -1 0 0 0; 0 -1 2 -1 0 0; 0 0 -1 2 -1 0; 0 0 0 -1 2 -1; 0 0 0 0 -1 2];

v0=rand(6,1);

v0 = v0/norm(v0);

for j=1:6

mu = 4*sin(j*pi/14)*sin(j*pi/14);

w = (T-mu*eye(6))\v0;

v = w/norm(w);

q = sqrt(2/7)*[sin(j*pi/7);sin(2*j*pi/7);sin(3*j*pi/7);sin(4*j*pi/7);sin(5*j*pi/7);sin(6*j*pi/7)];

norms = [norm(q-v) norm(q+v)]

end
