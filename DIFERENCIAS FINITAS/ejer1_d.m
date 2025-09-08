#ejercicio 1 - D)
clc;
clear;
#---------DATOS--------
k=1;
H=0.2;
c=1;
h=1/3;#Cant de puntos que se quiere
L1=0;
L2=1;
x=L1:h:L2;
n=length(x);
Q= 50;
Qinf=50;
K=zeros(n);
b=zeros(n,1);

#---------Diferencias centradas con nodo ficticio--------#
fila = [-1 (2+(h^2*c)/k) -1];
for i=2:(n-1)
  K(i,i-1:i+1)=fila;
end
K(1,1)=1;
K(n,n-1)=-2;
K(n,n)=(2 + (h^2*c)/k + (2*h*H/k));

b(1)=10;
b(2:n)=h*h*Q/k;
b(n)=b(n)+ (2*h*H*Qinf/k)

T_ap=K\b;

T_an=@(x) (-36.6897.*exp(-x) - 3.3103.*exp(x) + 50);


#-------GRAFICO
figure;
vec_x = linspace(L1,L2,n);
vec_y = T_an(vec_x);

plot(vec_x,T_ap,';T aproximado;');
hold on;
grid on;

plot(vec_x,vec_y,';T analitico;');
