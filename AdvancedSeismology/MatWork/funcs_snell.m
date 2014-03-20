function fs = funcs_snell(x,theta)
epsilon =0.2;
delta = 0.3;
Vp0 = 3;
Vs0 = 1.5;


f= 1- ((Vs0/Vp0)^2);
p1= 1+ (epsilon*((sind(x))^2)) - f/2;

p2= 4*(sind(x))^2*(1/f)*((2*delta*((cosd(x))^2)) - (epsilon*cosd(2*x)));

p3= 4* (epsilon^2) * ((sind(x))^4) *((1/f)^2);



V=p1 -((f/2)* sqrt(1 + p2 +p3));
V1 =Vp0*sqrt(V);

V1p0=3.0, V1s0=1.8, eps1=0.2, del1=0.1;
V2 = VSphase(theta,V1p0,V1s0,del1,eps1);


fs= (sind(x)/V1) - (sind(theta)/V2);




end
