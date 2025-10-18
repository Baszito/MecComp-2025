#ejercicio 1 - D)
#ACA TENER CONDICION ROBIN Y DIRICHLET
clc;
clear all
#---------DATOS--------
dx=4/10;
xnode=-4:dx:4;
k=3;
Q=100.*xnode;
H=2;
Qinf=10;
c=0;
#ARMADO MATRIZ K #BORDE ROBIN-DIRICHLET
#FUNCIONALIDAD FIJA

K=zeros(length(xnode),length(xnode));
b=zeros(length(xnode),1);
fila=[-k/dx^2 (2+c)*k/dx^2 -k/dx^2];
for i=2:length(xnode)-1
  K(i,i-1:i+1)=fila;
  b(i)=Q(i);
endfor
K(length(xnode),length(xnode)-1:length(xnode))=[-2*k/dx^2 2*k/dx^2 + 2*H/dx + c];
K(1,1)=1;

b(1)=20;
b(length(xnode))=Q(end)+(2*H/dx)*Qinf;

T=K\b
