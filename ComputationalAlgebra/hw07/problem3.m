format long e;
clear; clc;
A = zeros(64,64);
[n,m] = size(A);

for i=1:n
    A(i,i)=4;
end

for i=1:n-1
    A(i,i+1)=-1+1/100;
end

for i=2:n
    A(i,i-1)=-1-1/100;
end

e= ones(n,1);


k=30;
x_star= ones(n,1);

b=A*x_star;

Q=zeros(n,k+1);

Q(:,1)=b/norm(b);
H=zeros(k+1,k);
res = zeros(6,1);

for n=1:k
    v= A*Q(:,n);
    for i= 1:n
        H(i,n) = Q(:,i)'*v;
        v = v-H(i,n)*Q(:,i);
    end
    H(n+1,n)=norm(v);

    if H(n+1,n)== 0
        e1_exact=zeros(n,1);
        e1_exact(1,1)=1;
        y= H(1:n,1:n)\e1_exact;
        x_star_computed=norm(b)*Q(:,1:n)*y;
        break;
    end

    e1=zeros(n+1,1);
    e1(1,1)=1;
    [Q_h,R_h]=qr(H(1:n+1,1:n),0);
    y = R_h\(norm(b)*(Q_h)'*e1);
    x=Q(:,1:n)*y;
    Q(:,n+1)= v/(H(n+1,n));

    if rem(n,5)== 0
        norm(x_star-x)
    end
end