function out=Chua2(t,X)
% Values of parameters
alpha=9.0;
beta=100/7;
m0=-5/7;
m1=-8/7;

% Fuerza de conección
if t>10
c=0.0;
else 
    c=5.2;
end


% Matriz de coneciones
A=[-1, 1;
    1,-1];

% Nodo 1 
x1=X(1); 
x2=X(2); 
x3=X(3);

% Nodo 2
x4=X(4);
x5=X(5);
x6=X(6);

% Conección 
Conec=A*[x1;x4];

%Red de Chuas ecuaciones
% Del nodo 1
g1=m0*x1+(1/2)*(m1-m0)*(abs(x1+1)-abs(x1-1));
f1=alpha*(x2-x1-g1) + c*(-x1+x4);
f2=x1-x2+x3+ c*(-x2+x5); 
f3=-beta*x2+ c*(-x3+x6);

% Del nodo 2
g2=m0*x4+(1/2)*(m1-m0)*(abs(x4+1)-abs(x4-1));
f4=alpha*(x5-x4-g2) + c*(x1-x4);
f5=x4-x5+x6 +c*(x2-x5); 
f6=-beta*x5+  c*(x3-x6);
out=[f1;f2;f3;f4;f5;f6];

