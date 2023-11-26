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





%% Functions

function theta = DCM2EA313(C) 
    Cm = reshape(C,[3,3])' 

    theta(2) = acosd(Cm(3,3)) ;
    theta(1) = asind(Cm(3,1)/sind(theta(2))) ;
    theta(3) = acosd(Cm(2,3)/sind(theta(2))) ;
    
end 

function dx = DCMkinematics(t, C) 
    Cm = reshape(C,[3,3])' ;
   
    Omega = 20 ; % degree per second
    w = [Omega*sind(0.01*t) 0.01 Omega*cosd(0.02*t)] ; % omega-b/n vector
    wSkew = [ 0   -w(3) w(2) ;  %skew of omeg-b/n
              w(3) 0   -w(1) ;
             -w(2) w(1) 0    ] ;
    
    % Cbn dot 
    Cdotm = -wSkew*Cm ;
    Cdot = reshape(Cdotm',[1,9]);

    dx = Cdot' ;
end