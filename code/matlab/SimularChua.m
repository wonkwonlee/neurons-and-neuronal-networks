%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construction of the initial condition for the network

clear all
close all
clc

y0=[0.1,0.2,0.3];

[T,Y] = ode45(@Chua,[0, 67],y0); 

figure(1)
plot(T,Y)
figure(2)
plot3(Y(:,1),Y(:,2),Y(:,3))