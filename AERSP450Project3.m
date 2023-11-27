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

%% Part D
theta = zeros(length(Cbn),3) ;
for i = 1:length(Cbn)
    theta(i,:) = DCM2EA313(Cbn(i,:)) ;
end

f = figure ;
subplot(3,1,1)
plot(t,theta(:,1),'r')
title('3-1-3 Euler Angle 1 vs Time')
xlabel('Time (sec)')
ylabel('Euler Angle (degrees)')
subplot(3,1,2)
plot(t,theta(:,2),'g')
title('3-1-3 Euler Angle 2 vs Time')
xlabel('Time (sec)')
ylabel('Euler Angle (degrees)')
subplot(3,1,3)
plot(t,theta(:,3),'b')
title('3-1-3 Euler Angle 3 vs Time')
xlabel('Time (sec)')
ylabel('Euler Angle (degrees)')

%% Part E


%% Part F
beta1 = zeros(length(Cbn),4) ;
for i = 1:length(Cbn)
    beta1(i,:) = DCM2quat(Cbn(i,:)) ;
end

f = figure ;
subplot(4,1,1)
plot(t,beta1(:,1),'m')
title('Beta 0 vs Time')
xlabel('Time (sec)')
ylabel('Quaternion 4')
subplot(4,1,2)
plot(t,beta1(:,2),'r')
title('Beta 1 vs Time')
xlabel('Time (sec)')
ylabel('Quaternion 1')
subplot(4,1,3)
plot(t,beta1(:,3),'g')
title('Beta 2 vs Time')
xlabel('Time (sec)')
ylabel('Quaternion 2')
subplot(4,1,4)
plot(t,beta1(:,4),'b')
title('Beta 3 vs Time')
xlabel('Time (sec)')
ylabel('Quaternion 3')


%% Part G
beta2 = zeros(length(t),4) ;
quatErr = zeros(1,length(t)) ;
beta2(1,:) = [1 0 0 0] ;

for i = 2:length(t)
    Omega = 20 ;
    w = [Omega*sind(0.01*t(i-1)) 0.01 Omega*cosd(0.02*t(i-1))] ; % omega-b/n vector
    BwSkew = [0    -w(1) -w(2) -w(3) ;
              w(1)  0     w(3) -w(2) ;
              w(2) -w(3)  0     w(1) ;
              w(3)  w(2) -w(1)  0    ] ;
    
    phi = expm(BwSkew*dt/2) ;

    beta2(i,:) = phi * beta2(i-1,:)' ;

    quatErr(i) = norm(beta1(i,:)-beta2(i,:)) ;
end


%% Functions

function [C2,err] = discreteDCM(t,C0,dt,C1)  
    C2 = zeros(length(t),9) ;
    err = zeros(length(t),1) ;
    err(1) = 0 ;
    C2(1,:) = reshape(C0',[1,9]) ;
   
    for i = 2:length(t)
        CmPrev = reshape(C2(i-1,:),[3,3])' ;        
        Omega = 20 ; % degree per second
        w = [Omega*sind(0.01*t(i-1)) 0.01 Omega*cosd(0.02*t(i-1))] ; % omega-b/n vector
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

    theta(2) = acosd(Cm(3,3)) ;
    theta(1) = asind(Cm(3,1)/sind(theta(2))) ;
    theta(3) = acosd(Cm(2,3)/sind(theta(2))) ;
end 

function beta = DCM2quat(C) 
    Cm = reshape(C,[3,3])' ;

    beta(1) = sqrt((1 + Cm(1,1) + Cm(2,2) + Cm(3,3))) / 2 ;
    beta(2) = (Cm(2,3) - Cm(3,2)) / (4*beta(1)) ;
    beta(3) = (Cm(3,1) - Cm(1,3)) / (4*beta(1)) ;
    beta(4) = (Cm(1,2) - Cm(2,1)) / (4*beta(1)) ;
end 

function dx = DCMkinematics(t, C) 
    Cm = reshape(C,[3,3])' ;
   
    Omega = 20 ; % degree per second
    w = [Omega*sind(0.01*t) 0.01 Omega*cosd(0.02*t)] ; % omega-b/n vector
    wSkew = [ 0   -w(3) w(2) ;  %skew of omega-b/n
              w(3) 0   -w(1) ;
             -w(2) w(1) 0    ] ;
    
    % Cbn dot 
    Cdotm = -wSkew*Cm ;
    Cdot = reshape(Cdotm',[1,9]);

    dx = Cdot' ;
end

function dx = EAkinematics(t)
    
    dx = 0 ;
end