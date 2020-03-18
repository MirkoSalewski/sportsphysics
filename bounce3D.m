%3D bounce motion of an arbitrarily spinning ball
clear;

% gravity
g=9.81;      %[m/s^2]

%ball
m=0.43;       %mass [kg]
r=0.1;        %ball radius [m]
a=2/3;        %coefficient of I: 2/3 for hollow sphere, 2/5 for solid sphere
I=a*m*r^2;    %moment of inertia

%horizontal coefficients of restitution in x and z
%and vertical coefficient of restitution 
%ey is between 0 and 1
%ex is between -1 and 1. -1: no friction, 0: rolling, completely inelastic,
%1 completely elastic collision that reverse sign of the velocity of the
%contact point, between -1 and 0 is between gliding and rolling
ex=0.8;  
ey=.95;     
ez=ex;   

k=1; %bounce number

%start position
xb(k)=0;   
yb(k)=0;
zb(k)=0;

%velocity before the first bounce (before: b)
vxb(k)=1;
vyb(k)=-1;
vzb(k)=0;

%spin vector before first bounce
wxb(k)=0;
wyb(k)=0;
wzb(k)=10;

vvecb(:,k)=[vxb(k); vyb(k); vzb(k);...
    wxb(k); wyb(k); wzb(k)];

%vacuum solution of the parabola
hangtime0=2*vyb(k)/g;
height0=vyb(k)^2/(2*g);
rangex0=vxb(k)*hangtime0;
rangez0=vzb(k)*hangtime0;

tvac0 = linspace(0,hangtime0);
xvac0 = vxb(k)*tvac0;
yvac0 = vyb(k)*tvac0-g/2*tvac0.^2;
zvac0 = vzb(k)*tvac0;

%loop over bounces
kvec=1:20;

%velocity and spin after the kth bounce    
% vxa(k)=(1-a*ex)/(a+1)*vxb(k)-a*(ex+1)/(a+1)*r*wzb(k);
% vya(k)=-ey*vyb(k);
% vza(k)=(1-a*ez)/(a+1)*vzb(k)+a*(ez+1)/(a+1)*r*wxb(k);
% 
% wxa(k)=(a-ez)/(a+1)*wxb(k)+(ez+1)/(a+1)/r*vzb(k);
% wya(k)=wyb(k);
% wza(k)=(a-ex)/(a+1)*wzb(k)-(ex+1)/(a+1)/r*vxb(k);

%bounce matrix
BOUNCE=[(1-a*ex)/(a+1) 0 0 0 0 -a*(ex+1)/(a+1)*r;
    0 -ey 0 0 0 0;
    0 0 (1-a*ez)/(a+1) a*(ez+1)/(a+1)*r 0 0;
    0 0 (ez+1)/(a+1)/r (a-ez)/(a+1) 0 0;
    0 0 0 0 1 0;
    -(ex+1)/(a+1)/r 0 0 0 0 (a-ex)/(a+1)];

%air matrix, calculates conditions before bounce n+1 from after bounce n in vacuum
AIR=eye(6);
AIR(2,2)=-1;

for k=kvec
%calculate conditions after bounce from before bounce
vveca(:,k)=BOUNCE*vvecb(:,k);

%get velocities for easier coding
vx=vveca(1,k);
vy=vveca(2,k);
vz=vveca(3,k);

%calculate range in x and z and hangtime after the kth bounce
hangtime(k)=2*vy/g;
rangex(k)=vx*hangtime(k);
rangez(k)=vz*hangtime(k);

%calculate vacuum trajectoriy
tvac(:,k) = linspace(0,hangtime(k));
xvac(:,k) = xb(k)+vx*tvac(:,k);
yvac(:,k) = vy*tvac(:,k)-g/2*tvac(:,k).^2;
zvac(:,k) = zb(k)+vz*tvac(:,k);

%calculate conditions before (k+1)th bounce
%in vacuum only the sign of vy reverses
vvecb(:,k+1)=AIR*vveca(:,k);

%calculate new bounce point
xb(k+1)=xb(k)+rangex(k);
zb(k+1)=zb(k)+rangez(k);

end

%calculate the energy, the force integrals and verify the horizontal
%coeefficients of restitution
vxa=vveca(1,:);vya=vveca(2,:);vza=vveca(3,:);
vxb=vvecb(1,:);vyb=vvecb(2,:);vzb=vvecb(3,:);
wxa=vveca(4,:);wya=vveca(5,:);wza=vveca(6,:);
wxb=vvecb(4,:);wyb=vvecb(5,:);wzb=vvecb(6,:);


E=0.5*m*(vxa.^2+vya.^2+vza.^2+r^2*(wxa.^2+wya.^2+wza.^2));
FxdtF=m*(vxa-vxb(1:end-1));
FydtF=m*(vya-vyb(1:end-1));
FzdtF=m*(vza-vzb(1:end-1));

FxdtM=I*(wza-wzb(1:end-1))/r;
FzdtM=-I*(wxa-wxb(1:end-1))/r;

checkex=-(vxa+r*wza)./(vxb(1:end-1)+r*wzb(1:end-1));
checkez=-(vza-r*wxa)./(vzb(1:end-1)-r*wxb(1:end-1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figures
figure(10); clf; hold on; box on;
subplot(2,1,1); hold on; box on;
plot(xvac0,yvac0)
plot(xvac,yvac)
xlabel('x [m]')
ylabel('y [m]')
axis equal
title('x-y')
hold off
subplot(2,1,2); hold on; box on;
plot(zvac0,yvac0)
plot(zvac,yvac)
xlabel('z [m]')
ylabel('y [m]')
axis equal
title('z-y')
hold off

figure(20); clf; hold on; box on;
plot3(xvac0,yvac0,zvac0)
plot3(xvac,yvac,zvac)


figure(30);clf;hold on; box on;
subplot(2,3,1); box on;
plot(kvec,vxa,'o')
xlabel('bounce [-]')
ylabel('vx [m/s]')
subplot(2,3,2); box on;
plot(kvec,vza,'o')
xlabel('bounce [-]')
ylabel('vz [m/s]')
subplot(2,3,3); box on;
plot(kvec,wxa,'o')
xlabel('bounce [-]')
ylabel('wx [1/s]')
subplot(2,3,4); box on;
plot(kvec,wza,'o')
xlabel('bounce [-]')
ylabel('wz [1/s]')
subplot(2,3,5); box on;
plot(kvec,E,'o')
xlabel('bounce [-]')
ylabel('E/m [J/kg]')

figure(40);clf;
subplot(1,3,1); hold on; box on;
plot(kvec,FxdtF,'o')
plot(kvec,FxdtM,'x')
xlabel('bounce [-]')
ylabel('Fxdt [Ns]')
subplot(1,3,2); hold on;box on;
plot(kvec,FydtF,'o')
xlabel('bounce [-]')
ylabel('Fydt [Ns]')
subplot(1,3,3); hold on;box on;
plot(kvec,FzdtF,'o')
plot(kvec,FzdtM,'x')
xlabel('bounce [-]')
ylabel('Fzdt [Ns]')
