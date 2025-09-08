#ejercicio 1 - b)
clc;
clear;
#---------DATOS--------
k=1;
h=1/96;#Cant de puntos que se quiere
L1=1;
L2=5;
x=L1:h:L2;
n=length(x);
Q= 100*((x(1:n-1)-3).^2);
K=zeros(n);
b=zeros(n,1);

#---------Diferencias centradas con nodo ficticio--------#
fila = [1 -2 1];
for i=2:(n-1)
  K(i,i-1:i+1)=fila;
end
K(1,1)=-2;
K(1,2)=2;
K(n,n)=1;

b(1:n-1)=-(-Q./k)*(h^2);
b(1)=b(1)+4*h;
b(n)=0;

T_ap=K\b;

T_an=@(y) (-25.*y.^4 +300.*y.^3 -1350.*y.^2 +1906.*y + 2345)/3

#-------GRAFICO
figure;
vec_x = linspace(L1,L2,n);
vec_y = T_an(vec_x);

plot(vec_x,T_ap,';T aproximado;');
hold on;
grid on;

plot(vec_x,vec_y,';T analitico;');
