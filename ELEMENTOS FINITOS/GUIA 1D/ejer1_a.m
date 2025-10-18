#---------------GUIA 1D - FEM---------------#
#Ejercicio 1 - A
clc;
clear all;

#Datos:
pcp=0;
k=2;
c=0;
G=@(x) zeros(1,length(x))+100;

#Discretizacion
L1=0;
L2=1;
h=4;
xnode=L1:(L2-L1)/(h-1):L2;
G=G(xnode);
Fe=[  0.3333
   0.6667
   1.0000
   1.0000];
#Condiciones de borde :
T_L1=10;
T_L2=50;

#Armado de matriz K y vector f
K=zeros(h,h);
f=zeros(h,1);

fila=[-1 2 -1];

for i=2:h-1
  K(i,i-1:i+1)=k*(h-1)*fila;
  f(i)=100*(Fe(i));
endfor
K=[  84   -96    12     0;
   -96   276  -192    12;
    12  -192   360  -180;
     0    12  -180   168]
#imponemos condiciones de borde
K(1,1)=1;
f(1)=T_L1;

K(h,h)=1;
f(h)=T_L2;

#obtenemos la solucion aproximada
T_ap=K\f

#Solucion analitica, para verificar
T_an= @(y) -25.*y.^2 +65.*y +10;
T_an(xnode)

plot(xnode,T_ap,";4 Nodos;");
hold on;
grid on;

vec_x=linspace(L1,L2,h*4);
vec_y=T_an(vec_x);

plot(vec_x,vec_y,";Exacta(4*h nodos);");

