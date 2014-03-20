function Vp = VPphase(theta,Vp0,Vs0,delta,epsilon)
%Calculates the p-wave phase velocity
f= 1- ((Vs0/Vp0)^2);
p1= 1+ (epsilon*((sind(theta))^2)) - f/2;

p2= 4*(sind(theta))^2*(1/f)*((2*delta*((cosd(theta))^2)) - (epsilon*cosd(2*theta)));

p3= 4* (epsilon^2) * ((sind(theta))^4) *((1/f)^2);



Vp2=p1 +((f/2)* sqrt(1 + p2 +p3));
Vp =Vp0*sqrt(Vp2);
end

