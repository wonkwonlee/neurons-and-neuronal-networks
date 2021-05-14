function out=ChenNetODE(t,y0)
N=10;

%% States matrix
S=zeros(3,N);
for j=1:N
S(:,j)=y0(1+(j-1)*3:3*j);
end
%% The interconnection matrix

load MatrizWSN10 
%Acop=A;



%% Parameters
alpha=10.0;
beta=15.0;
gamma=0.0385;
a=-1.270;
b=-0.680;

c=0.0;
if t>80
c=8.05;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:N
x1=S(1,j);
x2=S(2,j);
x3=S(3,j);

%Coupling=diag([1,0,0])*((Acop(j,:)* S')');
if x1>1
    fx1=-b*x1-a+b;
end
if abs(x1)<=1
    fx1=-a*x1;
end
if x1<-1
    fx1=-b*x1+a-b;
end
 
yd((1+(j-1)*3:3*j))=[alpha*(x2-x1+fx1)+c*Acop(j,:)*S(1,:)';
                     x1-x2+x3;
                     -beta*x2-gamma*x3];
end
out=yd';