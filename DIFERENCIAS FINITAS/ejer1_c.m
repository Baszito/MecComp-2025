#ejercicio 1 - b)
clc;
clear all;
#---------DATOS--------
L1=1;
L2=5;
dx=1/64;
xnode=L1:dx:L2;
k=1;
G=100.*(xnode-3).^2;
q=2;


#armado matriz K y vector b de coeficientes
K=zeros(length(xnode),length(xnode));
b=zeros(length(xnode),1);
fila=[-k/dx^2 2*k/dx^2 -k/dx^2];
for i=2:length(xnode)-1
  K(i,i-1:i+1)=fila;
  b(i)=G(i);
endfor
K(length(xnode),length(xnode))=1;
K(1,1:2)=[2*k/dx^2 -2*k/dx^2];
b(1)+=2*q/dx;
b(end)=0;

T_ap=K\b;
T_an=@(y) (-25.*y.^4 +300.*y.^3 -1350.*y.^2 +1906.*y + 2345)/3;

#-------GRAFICO
T_an_evaluado = T_an(xnode);
figure(1);
plot(xnode,T_ap,'*-');
grid on;
hold on;
plot(xnode,T_an_evaluado,'*-');
