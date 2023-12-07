clc
clear
close all

totalT = 120 ; % 2 minutes in seconds
I = eye(3) ;   % 3x3 identity matrix
Cbn0 = I ;     % initial condition

%% Part A
dt = 0.01 ; % delta T = 0.01
init = reshape(Cbn0',[1,9]) ;
t = 0:dt:totalT ;
options = odeset('reltol',1e-12,'abstol',1e-12) ;
[t, Cbn] = ode45( @(t,C) DCMkinematics(t,C) , t , init, options) ;

orthCheck = zeros(1,length(Cbn)) ;
for i=1:length(Cbn) 
    Cm = reshape(Cbn(i,:),[3,3])' ;
    orthCheck(i) = norm(Cm*Cm' - I) ;
end

f = figure ;
subplot(1,1,1)
plot(t,orthCheck)
title('Norm of Orthogonality Constraint Violation vs Time')
xlabel('Time (sec)')
ylabel('Norm of Orthogonality Constraint Violation')

%% Part B
[CbnDiscrete1,err1] = discreteDCM(t,Cbn0,dt,Cbn) ;

f = figure ;
subplot(1,1,1)
plot(t,err1,'r')
hold on

%% Part C
dtC = 1 ;
tC = 0:dtC:totalT ;
CbnC = zeros(length(tC),9) ;
for i = 1:length(tC)
    CbnC(i,:) = Cbn((t==tC(i)),:) ;
end

[CbnDiscrete2,err2] = discreteDCM(tC,Cbn0,dtC,CbnC) ;

plot(tC,err2,'b')
title('Norm of Error vs Time')
xlabel('Time (sec)')
ylabel('Norm of Error')
legend('dt = 0.01','dt = 1')
hold off

%% Part D
theta1 = zeros(length(Cbn),3) ;
theta1(1,:) = [pi/2,0,-pi/2] ;          % known inital angles
nh = zeros(1,length(Cbn)) ;
n = 0 ;
for i = 2:length(Cbn)
    theta1(i,:) = DCM2EA313(Cbn(i,:)) ;
    if abs(theta1(i,1) - theta1(i-1,1)) > 1
        n = n + 1 ;
    end
    nh(i) = n ;
end
for i = 2:length(Cbn)
    theta1(i,1) = theta1(i,1) + 2*pi*nh(i) ;
end

f = figure ;
subplot(3,1,1)
plot(t,theta1(:,1),'r')
hold on
subplot(3,1,2)
plot(t,theta1(:,2),'r')
hold on
subplot(3,1,3)
plot(t,theta1(:,3),'r')
hold on

%% Part E
init = [pi/2,10^-10,-pi/2] ;
tEA = 0:dt:totalT ;
[tEA, theta2] = ode45( @(tEA, theta2) EAkinematics(tEA, theta2) , tEA , init, options) ;    % this method does not work because the inital condition is bascially at gimble lock since theta2 = 0

subplot(3,1,1)
plot(tEA,theta2(:,1),'b')
title('3-1-3 Euler Angle 1 vs Time')
xlabel('Time (sec)')
ylabel('Euler Angle (radians)')
legend('DCM2EA313','ODE45')
hold off
subplot(3,1,2)
plot(tEA,theta2(:,2),'b')
title('3-1-3 Euler Angle 2 vs Time')
xlabel('Time (sec)')
ylabel('Euler Angle (radians)')
legend('DCM2EA313','ODE45')
hold off
subplot(3,1,3)
plot(tEA,theta2(:,3),'b')
title('3-1-3 Euler Angle 3 vs Time')
xlabel('Time (sec)')
ylabel('Euler Angle (radians)')
legend('DCM2EA313','ODE45')
hold off

%% Part F
beta1 = zeros(length(Cbn),4) ;
quatCstr = zeros(length(Cbn),1) ;
for i = 1:length(Cbn)
    beta1(i,:) = DCM2quat(Cbn(i,:)) ;  
    quatCstr(i) = beta1(i,1)^2 + beta1(i,2)^2 + beta1(i,3)^2 + beta1(i,4)^2 ;
end

f = figure ;
subplot (1,1,1)
plot(t,quatCstr)
title('Quaternion Constraint vs Time')
xlabel('Time (sec)')
ylabel('Quaternion Constraint')
legend('dt = 0.01')

%% Part G
beta2 = zeros(length(t),4) ;
quatErr = zeros(length(t),1) ;
beta2(1,:) = [1 0 0 0] ;

for i = 2:length(t)
    Omega = 20 ;
    Omega = Omega * pi/180 ;
    w = [Omega*sin(0.01*t(i-1)) 0.01 Omega*cos(0.02*t(i-1))] ; % omega-b/n vector
    BwSkew = [0    -w(1) -w(2) -w(3) ;
              w(1)  0     w(3) -w(2) ;
              w(2) -w(3)  0     w(1) ;
              w(3)  w(2) -w(1)  0    ] ;
    
    phi = expm(BwSkew*dt/2) ;

    beta2(i,:) = phi * beta2(i-1,:)' ;
    
    beta2inv = [beta2(i,1) -beta2(i,2) -beta2(i,3) -beta2(i,4)] ;
    quatErr(i) = norm(qMult(beta1(i,:),beta2inv) - [1 0 0 0]') ;
    

    if quatErr(i) > 1                             % jank way of making the DCM results from Part F continuous
        beta1(i,:) = beta1(i,:) * -1 ;
        quatErr(i) = norm(qMult(beta1(i,:),beta2inv) - [1 0 0 0]') ;
    end    
end

f = figure ;
subplot(2,2,1)
plot(t,beta1(:,1),'r')
hold on
subplot(2,2,2)
plot(t,beta1(:,2),'r')
hold on
subplot(2,2,3)
plot(t,beta1(:,3),'r')
hold on
subplot(2,2,4)
plot(t,beta1(:,4),'r')
hold on
subplot(2,2,1)
plot(t,beta2(:,1),'b')
title('Beta 0 vs Time')
xlabel('Time (sec)')
ylabel('Quaternion 0')
legend('DCM2quat','Descretized')
hold off
subplot(2,2,2)
plot(t,beta2(:,2),'b')
title('Beta 1 vs Time')
xlabel('Time (sec)')
ylabel('Quaternion 1')
legend('DCM2quat','Descretized')
hold off
subplot(2,2,3)
plot(t,beta2(:,3),'b')
title('Beta 2 vs Time')
xlabel('Time (sec)')
ylabel('Quaternion 2')
legend('DCM2quat','Descretized')
hold off
subplot(2,2,4)
plot(t,beta2(:,4),'b')
title('Beta 3 vs Time')
xlabel('Time (sec)')
ylabel('Quaternion 3')
legend('DCM2quat','Descretized')
hold off

f = figure ;
subplot (1,1,1)
plot(t,quatErr)
title('Norm of Error vs Time')
xlabel('Time (sec)')
ylabel('Norm of Error')
legend('dt = 0.01')



%% Functions

function [C2,err] = discreteDCM(t,C0,dt,C1)  
    C2 = zeros(length(t),9) ;
    err = zeros(length(t),1) ;
    err(1) = 0 ;
    C2(1,:) = reshape(C0',[1,9]) ;
   
    for i = 2:length(t)
        CmPrev = reshape(C2(i-1,:),[3,3])' ;  
        Omega = 20 ; % degrees per second 
        Omega = Omega*pi/180 ;
        w = [Omega*sin(0.01*t(i-1)) 0.01 Omega*cos(0.02*t(i-1))] ; % omega-b/n vector
        wSkew = [ 0   -w(3) w(2) ;  %skew of omega-b/n
                  w(3) 0   -w(1) ;
                 -w(2) w(1) 0    ] ;
        
        phi = expm(-wSkew*dt) ;
        CmNew = phi*CmPrev ;
        C1m = reshape(C1(i,:),[3,3])' ; 
        err(i) = norm(C1m*CmNew' - eye(3)) ;
        C2(i,:) = reshape(CmNew',[1,9]) ;
    end
end

function theta = DCM2EA313(C) 
    Cm = reshape(C,[3,3])' ;
    
    theta(1) = atan2(Cm(3,1),-Cm(3,2)) ;
    theta(2) = acos(Cm(3,3)) ;
    theta(3) = atan2(Cm(1,3),Cm(2,3)) ;
end 

function beta = DCM2quat(C) 
    Cm = reshape(C,[3,3])' ;

    beta(1) = sqrt(Cm(1,1) + Cm(2,2) + Cm(3,3) + 1) / 2 ;
    beta(2) = (Cm(2,3) - Cm(3,2)) / (4*beta(1)) ; 
    beta(3) = (Cm(3,1) - Cm(1,3)) / (4*beta(1)) ; 
    beta(4) = (Cm(1,2) - Cm(2,1)) / (4*beta(1)) ;

end 

function res = qMult(a,b)
    res = [a(1)*b(1) - a(2)*b(2) - a(3)*b(3) - a(4)*b(4) ;
           a(1)*b(2) + a(2)*b(1) + a(3)*b(4) - a(4)*b(3) ;
           a(1)*b(3) - a(2)*b(4) + a(3)*b(1) + a(4)*b(2) ;
           a(1)*b(4) + a(2)*b(3) - a(3)*b(2) + a(4)*b(1) ] ;
end

function dx = DCMkinematics(t, C) 
    Cm = reshape(C,[3,3])' ;
   
    Omega = 20 ; % degree per second
    Omega = Omega*pi/180 ;
    w = [Omega*sin(0.01*t) 0.01 Omega*cos(0.02*t)]' ; % omega-b/n vector
    wSkew = [ 0   -w(3) w(2) ;  %skew of omega-b/n
              w(3) 0   -w(1) ;
             -w(2) w(1) 0    ] ;
    
    % Cbn dot 
    Cdotm = -wSkew*Cm ;
    Cdot = reshape(Cdotm',[1,9]);

    dx = Cdot' ;
end

function dx = EAkinematics(t, ang)
    theta1 = ang(1) ;
    theta2 = ang(2) ;
    theta3 = ang(3) ;

    Omega = 20 ; % degree per second
    Omega = Omega*pi/180 ;
    w = [Omega*sin(0.01*t) 0.01 Omega*cos(0.02*t)]' ; % omega-b/n vector
    
    B = 1/sin(theta2) * [sin(theta3) cos(theta3) 0 ; sin(theta2)*cos(theta3) -sin(theta2)*sin(theta3) 0 ; -cos(theta2)*sin(theta3) -cos(theta2)*cos(theta3) sin(theta2) ] ;
    
    thetaDot = B * w ;

    dx = thetaDot ;
end

% creates a dcm for an angle about axis 1
function r = dcm1axis(ang)
r = [1 0 0 ; 0 cos(ang) sin(ang) ; 0 -sin(ang) cos(ang)];
end

% creates a dcm for an angle about axis 3
function r = dcm3axis(ang)
r = [cos(ang) sin(ang) 0 ; -sin(ang) cos(ang) 0 ; 0 0 1];
end