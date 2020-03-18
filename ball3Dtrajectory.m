%This code calculates a 3D trajectory of an arbitrarily spinning ball
%subject to gravity, aerodynamic drag, and the Magnus force

clear; close all;

%physical parameters
g=9.81;       %[m/s^2]
rho=1.2;       % air density [kg/m^3]

%ball parameters
r=0.11;        %football radius [m]
Aproj=pi*r^2;  %projected or frontal area of the ball [m^2]
m=0.430;       %football kg

%default drag and lift coefficients 
cd=0.15;         %default ball drag coefficient
cl=1.7/pi;       %lift coefficient 

%initial position
x0=0;
y0=0;
z0=0;

%Roberto Carlos 1997 :13.5 rps, 38 m/s
%initial velocity
v0=38;                    %initial velocity in the (x,y) plane
theta0=45;                %initial angle to the horizontal
vz0=0;                    %initial velocity component in z direction
vtot0=sqrt(v0^2+vz0^2);   %initial total speed of the ball

%the spin vector is assumed to be constant during the flight
wx=0;
wy=2*pi*13.5;   
wz=0;

%magnitude of spin vector
w=sqrt(wx^2+wy^2+wz^2);


%non-dimensional parameters
A=v0^2/(r*g);                    %range parameter A
B=m*g/(0.5*cd*Aproj*rho*v0^2);   %ballistic parameter B
S=w*r/vtot0;                     %spin parameter S

%critical velocity for calculation of drag coefficient CD
vc=12.19; 
vs=1.309;

% here different loops are performed, which can be replaced by other loops
% such as %for m=[0.001 0.01 0.1 1 10], %for theta0=38:1:4, %for
% dto1=logspace(0,-3,4), or %for wy=[-10:2:10]*2*pi

%loop over different trajectories
%trajectory counter
kt=0;

%for v0=10:10:40
for v0=38 
kt=kt+1;
    
    %calculate velocity components
    vx0=cosd(theta0)*v0;
    vy0=sind(theta0)*v0;
    
    
    %characteristic parameters of the vacuum solution 
    hangtime_vac(kt)=2*vy0/g;
    height_vac(kt)=vy0^2/(2*g);
    range_vac(kt)=v0^2*sind(2*theta0)/g;
    
    %calculate analytic vacuum trajectory
    tvac(:,kt) = linspace(0,hangtime_vac(kt));
    xvac(:,kt) = vx0*tvac(:,kt);
    yvac(:,kt) = y0+vy0*tvac(:,kt)-g/2*tvac(:,kt).^2;
    
    
    %first-order numerical integration 
   
    %time step size (vary to check the convergence)
    dto1=0.001; %s
    
    %initial conditions
    xo1(1,kt)=x0;
    yo1(1,kt)=y0;
    zo1(1,kt)=z0;
    vxo1(1,kt)=vx0;
    vyo1(1,kt)=vy0;
    vzo1(1,kt)=vz0;
    to1(1,kt)=0;
    
    %step counter
    ko1=1;
    while min(yo1(:,kt))> -0.001  %a little less then zero to take out the initial point with z=0.
        
        %drag coefficient from v and S
        vo1(ko1,kt)=sqrt(vxo1(ko1,kt)^2+vyo1(ko1,kt)^2+vzo1(ko1,kt)^2);
        S=w*r/vo1(ko1,kt);
        if S>0.05 && vo1(ko1,kt)>vc
            cd=0.4127*S^0.3056;
            cdo1(ko1,kt)=0.4127*S^0.3056;
        else
            cd=0.155+0.346/(1+exp((vo1(ko1,kt)-vc)/vs));
            cdo1(ko1,kt)=0.155+0.346/(1+exp((vo1(ko1,kt)-vc)/vs));
        end
        
        %accelerations due to gravity, aerodyanmic drag and the Magnus force 
        axo1=-1/(2*m)*cd*Aproj*rho*vo1(ko1,kt)^2*vxo1(ko1,kt)/vo1(ko1,kt)...
            +cl/m*pi*r^3*rho*(wy*vzo1(ko1,kt)-wz*vyo1(ko1,kt));
        ayo1=-1/(2*m)*cd*Aproj*rho*vo1(ko1,kt)^2*vyo1(ko1,kt)/vo1(ko1,kt)-g...
            +cl/m*pi*r^3*rho*(wz*vxo1(ko1,kt)-wx*vzo1(ko1,kt));
        azo1=-1/(2*m)*cd*Aproj*rho*vo1(ko1,kt)^2*vzo1(ko1,kt)/vo1(ko1,kt)...
            +cl/m*pi*r^3*rho*(wx*vyo1(ko1,kt)-wy*vxo1(ko1,kt));
        
        %first-order integration of the acceleration gives the velocity,
        %and another first-order integration gives the position
        vxo1(ko1+1,kt)=vxo1(ko1,kt)+axo1*dto1;
        vyo1(ko1+1,kt)=vyo1(ko1,kt)+ayo1*dto1;
        vzo1(ko1+1,kt)=vzo1(ko1,kt)+azo1*dto1;
        xo1(ko1+1,kt)=xo1(ko1,kt)+vxo1(ko1,kt)*dto1;
        yo1(ko1+1,kt)=yo1(ko1,kt)+vyo1(ko1,kt)*dto1;
        zo1(ko1+1,kt)=zo1(ko1,kt)+vzo1(ko1,kt)*dto1;
        
        %advance the time
        to1(ko1+1,kt)=to1(ko1,kt)+dto1;
        
        ko1=ko1+1;
    end
    
    %numerically calculated characteristic trajectory parameters
    rangexo1(kt) = max(xo1(:,kt));
    rangezo1(kt) = max(zo1(:,kt));
    heighto1(kt) = max(yo1(:,kt));
    hangtimeo1(kt)=max(to1(:,kt));
    
    
    %this section uses the Runge-Kutta scheme implemented in Matlab for the trajectory calculation    
    tspan=[0 1.1*hangtimeo1(kt)];
    eq1=@(t,x)[x(2);1/m*(-1/2*cdfunc(x(2),x(4),x(6),w,r,vc,vs)*Aproj*rho*(x(2)^2+x(4)^2+x(6)^2)*x(2)/(sqrt(x(2)^2+x(4)^2+x(6)^2))...
        +cl*pi*r^3*rho*(wy*x(6)-wz*x(4)));...
        x(4);1/m*(-m*g-1/2*cdfunc(x(2),x(4),x(6),w,r,vc,vs)*Aproj*rho*(x(2)^2+x(4)^2+x(6)^2)*x(4)/(sqrt(x(2)^2+x(4)^2+x(6)^2))...
        +cl*pi*r^3*rho*(wz*x(2)-wx*x(6)));...
        x(6);1/m*(-1/2*cdfunc(x(2),x(4),x(6),w,r,vc,vs)*Aproj*rho*(x(2)^2+x(4)^2+x(6)^2)*x(6)/(sqrt(x(2)^2+x(4)^2+x(6)^2))...
        +cl*pi*r^3*rho*(wx*x(4)-wy*x(2)))];
    
    sol1=ode45(eq1,tspan,[x0 vx0 y0 vy0 z0 vz0]);
    trk45 = linspace(0,1.1*hangtimeo1(kt));
    xrk45(:,kt) = deval(sol1,trk45,1);
    vxrk45(:,kt)= deval(sol1,trk45,2);
    yrk45(:,kt) = deval(sol1,trk45,3);
    vyrk45(:,kt)= deval(sol1,trk45,4);
    zrk45(:,kt) = deval(sol1,trk45,5);
    vzrk45(:,kt)= deval(sol1,trk45,6);
    
    xrk45(find(yrk45(:,kt)<-eps),kt)=nan;
    vxrk45(find(yrk45(:,kt)<-eps),kt)=nan;
    vyrk45(find(yrk45(:,kt)<-eps),kt)=nan;
    zrk45(find(yrk45(:,kt)<-eps),kt)=nan;
    vzrk45(find(yrk45(:,kt)<-eps),kt)=nan;
    yrk45(find(yrk45(:,kt)<-eps),kt)=nan;
end

%the num1 solution has different lengths for different trajectories,
%we replace last zeros with nan to prevent plotting

[npoints,ntrajectories]=size(xo1);
for kp=2:npoints
    for kt=1:ntrajectories
        if xo1(kp,kt)<eps
            xo1(kp,kt)=nan;
            yo1(kp,kt)=nan;
            zo1(kp,kt)=nan;
        end
    end
end



%compare vacuum, the order 1, and the Runge-Kutta solutions
figure(10); clf; 
subplot(2,1,1); hold on; box on;
plot(xrk45,yrk45)
plot(xo1,yo1)
plot(xvac,yvac)
xlabel('x [m]')
ylabel('y [m]')
xlim([0 1.1*max(max(max(xvac)),max(max(xo1)))])
ylim([0 1.1*max(max(max(yvac)),max(max(yo1)))])
title('x-y')
hold off
subplot(2,1,2); hold on; box on;
plot(xrk45,zrk45)
plot(xo1,zo1)
xlabel('x [m]')
ylabel('z [m]')
axis equal
title('x-z')
hold off


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %drag coefficient function for Runge-Kutta integration
function cd=cdfunc(vx,vy,vz,w,radius,vc,vs)
v=sqrt(vx^2+vy^2+vz^2);
S=w*radius/v;
        if S>0.05 && v>vc
            cd=0.4127*S^0.3056;
        else
            cd=0.155+0.346/(1+exp((v-vc)/vs));
        end

end
