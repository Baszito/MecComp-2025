#ejercicio 1 - A)
clc;
clear;
#---------DATOS--------
k=2;
c=0;
Q=100;
h=1/3;
L1=0;
L2=1;
x=L1:h:L2;
n=length(x)
K=zeros(n,n);
b=zeros(n,1);

fila = [1 -2 1];
for i=2:(n-1)
  K(i,i-1:i+1)=fila;
end
K(1,1)=1;
K(n,n)=1;

b(1)=10;
b(2:n-1)=-(Q./k)*(h^2);
b(n)=50;

T_ap=K\b

T_an= @(y) -25.*y.^2 +65.*y +10;
T_an(x)

vec_x = linspace(L1,L2,n);
vec_y = T_an(vec_x)

plot(vec_x,T_ap);
hold on;
grid on;

plot(vec_x,vec_y);
