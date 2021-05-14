function out=Chua(t,y)
alpha=10.0;
beta=15.0;
gamma=0.0385;
a=-1.270;
b=-0.680;

if t>25
T=125.0;
else
    T=0;
end


x1=y(1);
x2=y(2);
x3=y(3);

if x1>1
    fx1=-b*x1-a+b;
end
if abs(x1)<=1
    fx1=-a*x1;
end
if x1<-1
    fx1=-b*x1+a-b;
end

out=[alpha*(x2-x1+fx1)-T*x1;
    x1-x2+x3;
    -beta*x2-gamma*x3];