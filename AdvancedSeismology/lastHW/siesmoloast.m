clear all;close all;
%
% Nonhyperbolic (long spread) reflection moveout of P waves in layered VTI 
% media.

% Equation parameters.

% Interval parameters.

%  depth     Vp0        epsilon      delta
%  (km)     (km/s)
%  1.480    1.941        0.230       0.011
%  1.980    2.168        0.247       0.010
%  2.244    2.685        0.264       0.016

% Data provided in the home assignment.

load Traveltimes.dat;
%x - offset
global x 
x = Traveltimes(:,1);

%-traveltime
global torig
torig = Traveltimes(:,2);

% Depths (Km).  
depth1 = 1.480;
depth2 = 1.980;
depth3 = 2.244;

% P-wave vertical velocities (Km/s).
Vp01 = 1.941;
Vp02 = 2.168;
Vp03 = 2.685;

% Epsilons.
epsilon1 = 0.230;
epsilon2 = 0.247;
epsilon3 = 0.264;

% Deltas.
delta1 = 0.011;
delta2 = 0.010;
delta3 = 0.016;

% Etas.
eta1 = (epsilon1-delta1)/(1+(2*delta1));
eta2 = (epsilon2-delta2)/(1+(2*delta2));
eta3 = (epsilon3-delta3)/(1+(2*delta3));

% Zero-offset traveltimes (s).
t01 = 2*(depth1/Vp01);
t02 = 2*((depth2-depth1)/Vp02);
t03 = 2*((depth3-depth2)/Vp03);

% NMO velocities (Km/s). Vnmo(0)
Vnmo1 = Vp01*sqrt(1+(2*delta1));
Vnmo2 = Vp02*sqrt(1+(2*delta2));
Vnmo3 = Vp03*sqrt(1+(2*delta3));

% Effective parameters

global t0
t0 = (t01+t02+t03);
global Vnmo
%using dix type equation averaging
Vnmo = sqrt((1/t0)*(((Vnmo1^2)*t01)+((Vnmo2^2)*t02)+((Vnmo3^2)*t03)));
global eta
eta = (1/8)*((1/((Vnmo^4)*t0))*(((Vnmo1^4)*(1+(8*eta1))*t01)+((Vnmo2^4)*(1+(8*eta2))*t02)+((Vnmo3^4)*(1+(8*eta3))*t03))-1);
global Vhor
Vhor = Vnmo*sqrt(1+(2*eta));

% Nonhyperbolic moveout equation (effective quantities).

 C = 1.1;

t9 = sqrt(t0+((x.^2)/(Vnmo^2))-(((Vhor^2)-(Vnmo^2)).*(x.^4)./((Vnmo^2).*(((t0^2).*(Vnmo^4))+(C.*(Vhor^2).*(x.^2))))));

%Least Squares
%L=length(x);
k=100
C=zeros(k,1);
d1=0;
for j=1:k
    C(j)=0.8 + ((j-1)/200);
    d(j)=0;
 
    for i=1:length(x)
        p1(i,j)=(t0^2)+((x(i).^2)/(Vnmo^2));
        p2(i,j)=((Vhor^2)-(Vnmo^2))*(x(i))^4;
        p3(i,j)=(((t0^2).*(Vnmo^4)) + (C(j).*(Vhor^2).*(x(i).^2)))*(Vnmo^2);
        t(i,j)= sqrt(p1(i,j)- (p2(i,j)/p3(i,j)));
    
        d(j) = d(j)+[ t(i,j)-torig(i)]^2;
        mitr(i,j)= [ t(i,j)-torig(i)]^2;
    end
    
    if j==1
       d1=d(j);
    end
    if ( (d(j) <d1) )
        d1=d(j);
        minj=j-1 
    else j=1;
    end

end

M=zeros(1,k);
for j=1:k
    for i=1:length(x)
        M(j)=M(j)+mitr(i,j);
    end
end


%figure(1);
%t1=t(:,1)
%plot(t1(40:50));
figure(1);
plot(t(:,40),'o');hold on; plot(t(:,minj),'r');plot(torig,'+');ylim([2.7 3.1]);legend('c=1','c=1.155','exact traveltime','location','Best')
 xlabel('offset');title('Travel-times for various values of C')
ylabel('t(s)');

%mp=1:0.005:1.2.995;

figure(2);
plot(C,M); xlabel('C');title('Least squares error for various values of C from 0.8 - 1.3')
ylabel('least squares error');%set(gca,'XTick',[1:0.005:1.3])


d3=zeros(k,k);

%SECOND part : finding best fit eta
for j=1:k
    Vhor_o(j)= Vhor*((0.01*j) +0.49);
    for z=1:k
        Vnmo_o(z)= Vnmo*((0.01*z) +0.49);
        %d3(j,z)=0;
    
        for i=1:length(x)

            p1(i,j,z)=(t0^2)+((x(i).^2)/(Vnmo_o(z)^2));
            p2(i,j,z)=((Vhor_o(j)^2)-(Vnmo_o(z)^2))*(x(i))^4;
            p3(i,j,z)=(((t0^2).*(Vnmo_o(z)^4)) + (1.*(Vhor_o(j)^2).*(x(i).^2)))*(Vnmo_o(z)^2);
            t3(i,j,z)= sqrt(p1(i,j,z)- (p2(i,j,z)/p3(i,j,z)));
    
            d3(j,z) = d3(j,z)+[ t3(i,j,z)-torig(i)]^2;
            mitr(i,j,z)= [ t3(i,j,z)-torig(i)]^2;
        end
    
    if (j==1 && z==1)
       d1=d3(j,z);
    end
    if ( (d3(j,z) <d1) && j~=1)
        d1=d3(j,z);
        vn=Vnmo_o(z);
        vh=Vhor_o(j);
        minj3=j
        minz3=z
    end
    end
end

figure(3);surf(d3)
figure(4);
plot(t(:,40),'o');hold on; plot(t(:,minj),'r');plot(torig,'+');plot(t3(:,50,50),'-');ylim([2.7 3.1]);legend('c=1','c=1.155','exact traveltime','Vnmo and Vhor varying','location','Best')
 xlabel('offset');title('Travel-times')
ylabel('t(s)');
%ete_right=0.5*(())
% Least squares function.

%    Cinit = 1;  % Make a starting guess at the solution (C = 1).   
%    options = optimset('Display','iter');   % Option to display output.
%    [C] = lsqnonlin(@objfun,Cinit,options);  % Call solver.