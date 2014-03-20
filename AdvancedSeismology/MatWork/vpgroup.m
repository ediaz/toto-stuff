function [VGp, Psip] = vpgroup( j,Vp0,Vs0,delta,epsilon)
%vgroup and angle


%Vp=VPphase(theta,Vp0,Vs0,delta,epsilon)

  %  if j==2 || j==89
   %     dVp= (VPphase(j+1,Vp0,Vs0,delta,epsilon)-VPphase(j-1,Vp0,Vs0,delta,epsilon))/2;
   % else
   %dVp= (-VPphase(j+2,Vp0,Vs0,delta,epsilon) + (8*VPphase(j+1,Vp0,Vs0,delta,epsilon)) - (8*VPphase(j-1,Vp0,Vs0,delta,epsilon)) + VPphase(j-2,Vp0,Vs0,delta,epsilon))/12;
    
   f= 1- ((Vs0/Vp0)^2);
   Vp= VPphase(j,Vp0,Vs0,delta,epsilon);
    
   dVp=(Vp0^2*0.5*epsilon*sind(2.*j))./Vp + (Vp0^2./Vp).*((0.5*sind(2.*j).*(2*delta*(cosd(j).*cosd(j)) - epsilon*cosd(2.*j)) + sind(j).*sind(j).*(epsilon*sind(2.*j) - delta* sind(2.*j)) + (2/f)*epsilon^2.*(sind(j).^3.*cosd(j)))./(sqrt((1. + (4/f).*(sind(j).^2).*(2*delta*cosd(j).^2 - epsilon*cosd(2.*j)) + ((2*epsilon/f)^2).*(sind(j).^4))))); 
  
  
   
    
    
    
    VGp= Vp* sqrt(1 + ((dVp/Vp)^2));
    
    Psip= atand((tand(j) + (dVp/Vp)) / (1- (tand(j)*(dVp/Vp))));
    


end

