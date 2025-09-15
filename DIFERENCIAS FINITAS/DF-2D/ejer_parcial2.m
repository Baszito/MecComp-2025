#Otro script
clc;
clear all;
#---------DATOS--------
dx=4/10;
xnode=-4:dx:4;
k=3;
c=0;
G=100.*(xnode);
h=2;
Tinf=10;

#armado matriz K y vector b de coeficientes
K=zeros(length(xnode),length(xnode));
b=zeros(length(xnode),1);
fila=[-k/dx^2 2*k/dx^2 -k/dx^2];
for i=2:length(xnode)-1
  K(i,i-1:i+1)=fila;
  b(i)=G(i);
endfor

#CONDICION DIRICHLET BORDE IZQUIEDO
K(1,1)=1;
b(1)=20;


#CONDICION ROBIN BORDE DERECHO
K(length(xnode),length(xnode)-1:length(xnode))=[-2*k/dx^2 (2+(2*dx*h/k))*k/dx^2];
b(end)=G(xnode(end)) - 2*h*Tinf/dx;

T_ap=K\b;

