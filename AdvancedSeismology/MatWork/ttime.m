%Siesmo three travel times of reflected P and S waves
clear all; close all;

eps=0.2;
del= 0.1;
Vp0=3.0;
Vs0=1.8;
h1=0.5;
h2=1.0;
teta1=1:1:90;

for i=1:90;
    theta=i-1;
    
    Vp(i)= VPphase(theta,Vp0, Vs0, del, eps);
    Vs(i)= VSphase(theta,Vp0, Vs0, del, eps);
    
    %Vx(i)= VPphase(x,Vp0, Vs0, del, eps);
    %eqn:=x=asind((sind(theta)*VPphase(x,3.0,1.5,0.3,0.2))/Vp(i))
    
    %angI=theta;
    
    options=optimset('MaxFunEvals',500,'MaxIter',50,'TolFun',1e-12,'Display','off');
    [angPsnell,FVAL,EXITFLAG] = fzero(@(x) funcp_snell(x,theta),theta,options);
    
    angP_t(i)=abs(angPsnell);
   % x=fzero(@(x) (sind(theta)/Vp(i)-(sind(x)/VPphase(x,Vp0, Vs0, del, eps))),x)
    
    [VGp,Psip]=vpgroup(theta,Vp0,Vs0,del,eps);
    Vg1(i)=VGp;
    PsiP1(i)=Psip;
    
    [VGp2,Psip2]=vpgroup(angP_t(i),Vp0,1.5,0.3,-0.2);   
    Vg2(i)=VGp2;
    PsiP2(i)=Psip2;
    
    %layer 1
    off1(i)=0.5*tand(Psip);
    t1(i)  = (0.5/cosd(Psip))/VGp;
    
    %layer 2
    off2(i)=1*tand(Psip2);
    t2(i)  = (1/cosd(Psip2))/VGp2;    
    
    offp= off1(i)+off2(i)
    if abs(offp) <= 1
       offset(i)=2* (off1(i)+off2(i)); 
       traveltime(i)= 2*(t1(i)+ t2(i));
       off_neg(i)=-offset(i);
    end
    
    %offset(i)=2* (off1(i)+off2(i));
    %traveltime(i)= 2*(t1(i)+ t2(i));
end

figure;

%VGp=vpgroup(0,Vp0,Vs0,del,eps);
plot(offset,traveltime);

hold on;
%VGp=vpgroup(0,Vp0,Vs0,del,eps);
plot(off_neg,traveltime);hold on;
title('Pwave, two-way traveltime for epsilon=-0.2');
 set(gca,'YDir','Reverse')
 xlabel('Offset(km)')
ylabel('t(s)');
 
for i=1:90;
    theta=i-1;
    %S-waves
    options=optimset('MaxFunEvals',500,'MaxIter',50,'TolFun',1e-12,'Display','off');
    [angSsnell,FVAL,EXITFLAG] = fzero(@(x1) funcs_snell(x1,theta),theta,options);
    
    angS_t(i)=abs(angSsnell);
   % x=fzero(@(x) (sind(theta)/Vp(i)-(sind(x)/VPphase(x,Vp0, Vs0, del, eps))),x)
    
    [VGs,Psis]=vsgroup(theta,Vp0,Vs0,del,eps);
    Vgs1(i)=VGs;
    PsiS1(i)=Psis;
    
    [VGs2,Psis2]=vsgroup(angS_t(i),Vp0,1.5,0.3,0.2);   
    Vgs2(i)=VGs2;
    PsiS2(i)=Psis2;
    
    Slowness(i)=1/VSphase(angS_t(i),Vp0,1.5,0.3,0.2);
    Scos(i)=Slowness(i)*cosd(angS_t(i));
    Ssin(i)=Slowness(i)*sind(angS_t(i));
    
    
    %layer 1
    off1s(i)=0.5*tand(Psis);
    %t1s(i)  = (0.5/cos(Psis))/VGs;
    
    %layer 2
    off2s(i)=1*tand(Psis2);
    %t2s(i)  = (1/cos(Psis2))/VGs2;    
    
    
    offsets= (off1s(i)+off2s(i));
    if abs(offsets) <= 1 
        cmp1(i)=off1s(i);
        cmp2(i)=off2s(i);
   %     traveltimes(i)= 2*(t1s(i)+ t2s(i));
        offset_s(i)=2*(cmp1(i)+cmp2(i));
        travel_time_s(i)= 2*sqrt((cmp1(i)^2)+0.25)/Vgs1(i) + 2*sqrt((cmp2(i)^2)+1)/Vgs2(i);
        off_negs(i)=-offset_s(i);
    end
end


%figure;


figure;
%VGp=vpgroup(0,Vp0,Vs0,del,eps);
hold on;
plot(offset_s,travel_time_s);
plot(off_negs,travel_time_s); hold on;


title('Shear wave, two-way traveltime for epsilon=0.2');
 set(gca,'YDir','Reverse')
 xlabel('Offset(km)')
ylabel('t(s)');
 
figure;
polar(Ssin,Scos);

figure;
hold on;
for i=1:90
if (isnan(cmp1(i)) == 0)
        
     plot([0 cmp1(i)], [0 0.5]);
     plot([cmp1(i) 2*cmp1(i)], [0.5 0]);
        
end   
    
if (isnan(cmp1(i)) == 0 && isnan(cmp2(i)) == 0)
        
     plot([0 cmp1(i)], [0 0.5]);
     plot([cmp1(i) (cmp1(i)+cmp2(i))], [0.5 1.5])
     plot([(cmp1(i)+cmp2(i)) (cmp1(i)+2*cmp2(i))], [1.5 0.5])
     plot([(cmp1(i)+2*cmp2(i)) 2*(cmp1(i)+cmp2(i))], [0.5 0])
        
end
set(gca,'YDir','Reverse')
set(gca,'XLim',[0 2])
title('Shear wave, two-way traveltime for epsilon=-0.2');

 xlabel('Offset(km)')
ylabel('t(s)');

end

for i=1:1:90

    theta1 = PsiS1(i);
    theta2 = PsiS2(i);
    
    %find the corresponding cmps
    
    
    if (isinf(tand(theta1)) == 0)
    cmp1 = 0.5*tan(theta1);
    end
    
    if (isinf(tand(theta2)) == 0)
        cmp2 = tand(theta2);
    end
    
     if ((cmp1 + cmp2) <= 1)
        
        C(4,i) = cmp1;
        C(5,i) = cmp2;
        
     end
    
end


for i=1:1:90

    theta1 = PsiS1(i);
    
    
    %find the corresponding cmps
    
    
    if (isinf(tand(theta1)) == 0)
    cmp1 = 0.5*tand(theta1);
    end
    
    
     if (cmp1  <= 1)
        
        C(6,i) = cmp1;
        
        
     end
    
end

for i =1:1:90
    

    
    t4 = 2*sqrt(C(4,i)^2 + 0.25)/Vgs1(i) + 2*sqrt(C(5,i)^2 + 1)/Vgs2(i);
    d4 = 2*C(4,i) + 2*C(5,i);
    
       
    ts2(1,i) = d4;
    ts2(2,i) = t4;
    
end



figure;
plot(ts2(1,:),ts2(2,:));



