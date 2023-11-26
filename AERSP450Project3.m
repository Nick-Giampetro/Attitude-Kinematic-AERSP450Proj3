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




%% Functions

function dx = DCMkinematics(t, C) 
    Cm = reshape(C,[3,3])' ;
   
    Omega = 20 ; % degree per second
    w = [Omega*sind(0.01*t) 0.01 Omega*cosd(0.02*t)] ;
    wSkew = [ 0   -w(3) w(2) ;
              w(3) 0   -w(1) ;
             -w(2) w(1) 0    ] ;
    
    Cdotm = -wSkew*Cm ;
    Cdot = reshape(Cdotm',[1,9]);

    dx = Cdot' ;
end