%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construction of the initial condition for the network

clear all
close all
clc

N=10; %Number of nodes

y01=[0.1,0.2,0.3]; %Initial conditions Lorenz Equation 

% For the entire network the initial condition is given by
y0=zeros(1,N*3);
y0(1:3)=y01; % The first node is the reference 
             % and has the original initial conditions

rand('state',0);
for j=2:N
    y0(1+(j-1)*3:3*j)=y01*rand*2;
end
% y0 is the initial condition for all the nodes in a row-vector

[T,Y] = ode23(@ChuaNetODE,[0, 255],y0); % Solution of the ODE LorenzNet 
plot(T,Y,'k')
