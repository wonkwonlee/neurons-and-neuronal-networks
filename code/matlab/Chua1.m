function out=Chua1(t,X)
% Values of parameters
alpha=9.0;
beta=100/7;
m0=-5/7;
m1=-8/7;

x1=X(1); 
x2=X(2); 
x3=X(3);

%Chua equation
g=m0*x1+(1/2)*(m1-m0)*(abs(x1+1)-abs(x1-1));
f1=alpha*(x2-x1-g);
f2=x1-x2+x3; 
f3=-beta*x2;

out=[f1;f2;f3];

