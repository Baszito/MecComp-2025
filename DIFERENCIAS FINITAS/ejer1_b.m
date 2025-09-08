#ejercicio 1 - b)
clc;
clear;
#---------DATOS--------
k=1;
c=1;
Q=0;
h=1/3;
L1=0;
L2=2;
x=L1:h:L2;
n=length(x)
K=zeros(n,n);
b=zeros(n,1);

#---------Diferencias no centradas-hacia atras---------#
fila = [-1 (2+(h^2*c)/k) -1];
for i=2:(n-1)
  K(i,i-1:i+1)=fila;
end
K(1,1)=1;
K(n,n-1)=-1;
K(n,n)=1;

b(1)=100;
b(2:n-1)=-(Q./k)*(h^2);
b(n)=0;

T_ap=K\b;

T_an= @(y) ((100*(e.^(-y))).*((e.^(2.*y))+e^4))./(1+e^4);
T_an(x);

#-------GRAFICO
figure;
vec_x = linspace(L1,L2,n);
vec_y = T_an(vec_x);

plot(vec_x,T_ap);
hold on;
grid on;

plot(vec_x,vec_y);

#---------Diferencias centradas con nodo ficticio--------#
fila = [-1 (2+(h^2*c)/k) -1];
for i=2:(n-1)
  K(i,i-1:i+1)=fila;
end
K(1,1)=1;
K(n,n-1)=-2;
K(n,n)=(2+(h^2*c)/k);

b(1)=100;
b(2:n-1)=-(Q./k)*(h^2);
b(n)=0;

T_ap=K\b;

T_an= @(y) ((100*(e.^(-y))).*((e.^(2.*y))+e^4))./(1+e^4);
T_an(x);

#-------GRAFICO
figure;
vec_x = linspace(L1,L2,n);
vec_y = T_an(vec_x);

plot(vec_x,T_ap);
hold on;
grid on;

plot(vec_x,vec_y);

