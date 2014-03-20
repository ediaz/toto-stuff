function [VGs, Psis] = vsgroup( j,Vp0,Vs0,delta,epsilon)
%vgroup and angle


%Vp=VPphase(theta,Vp0,Vs0,delta,epsilon)

  %  if j==2 || j==89
  %      dVs= (VSphase(j+1,Vp0,Vs0,delta,epsilon)-VSphase(j-1,Vp0,Vs0,delta,epsilon))/2;
  %  else
      %  dVs= (-VSphase(j+2,Vp0,Vs0,delta,epsilon) + (8*VSphase(j+1,Vp0,Vs0,delta,epsilon)) - (8*VSphase(j-1,Vp0,Vs0,delta,epsilon)) + VSphase(j-2,Vp0,Vs0,delta,epsilon))/12;
    f= 1- ((Vs0/Vp0)^2);
      Vs=VSphase(j,Vp0,Vs0,delta,epsilon)      
   dVs= (Vp0^2*0.5*epsilon*sind(2.*j))./Vs - (Vp0^2./Vs).*((0.5*sind(2.*j).*(2*delta*(cosd(j).*cosd(j)) - epsilon*cosd(2.*j)) + sind(j).*sind(j).*(epsilon*sind(2.*j) - delta* sind(2.*j)) + (2/f)*epsilon^2.*(sind(j).^3.*cosd(j)))./(sqrt((1. + (4/f).*(sind(j).^2).*(2*delta*cosd(j).^2 - epsilon*cosd(2.*j)) + ((2*epsilon/f)^2).*(sind(j).^4)))));
   
  
    
   
    
    VGs= Vs* sqrt(1 + ((dVs/Vs)^2));
    
    Psis= atand((tand(j) + (dVs/Vs)) / (1- (tand(j)*(dVs/Vs))));
    
    
    
end

