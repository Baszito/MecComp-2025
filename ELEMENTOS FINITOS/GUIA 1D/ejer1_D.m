#---------------GUIA 1D - FEM---------------#
#Ejercicio 1 - A
clc;
clear all;

#Datos:
pcp=0;
k=1;
c=1;
G=@(x) zeros(length(x),1)+50;
H=0.2;

#Discretizacion
L1=0;
L2=1;
cant_nodos=4;
h=(L2-L1)/(cant_nodos-1);
xnode=L1:h:L2;
G=G(xnode);

#Condiciones de borde :
T_L1=10;
qinf_L2=50;

#Armado de matriz Ke

#Aporte Difusivo :
De=zeros(2,2);
De(1,1)=k/h;
De(1,2)=-k/h;
De(2,1)=-k/h;
De(2,2)=k/h;

#Aporte Reactivo
Ce=zeros(2,2);
Ce(1,1)=2*c*h/6;
Ce(1,2)=c*h/6;
Ce(2,1)=c*h/6;
Ce(2,2)=2*c*h/6;



#Armado de matriz K
K=zeros(cant_nodos,cant_nodos);
F=zeros(cant_nodos,1);
for i=1:cant_nodos-1
  K(i,i)+=De(1,1)+Ce(1,1);
  K(i+1,i)+=De(2,1)+Ce(2,1);
  K(i,i+1)+=De(1,2)+Ce(1,2);
  K(i+1,i+1)+=De(2,2)+Ce(2,2);

  #Armado vector F
  F(i)+=G(i)*h/2;
  F(i+1)+=G(i+1)*h/2;
endfor

#imponemos condiciones de borde

#Condicion Dirichlet :
K(1,:)=zeros(1,cant_nodos);
K(1,1)=1;
F(1)=T_L1;

#Condicion Neumann :
F(cant_nodos)+=qinf_L2*H;
#Aporte Robin
K(cant_nodos,cant_nodos)+=H;

F
#obtenemos la solucion aproximada
T_ap=K\F

#Solucion analitica, para verificar
T_an= @(y) -(36.6897.*e.^(-y)) -(3.3103.*e.^y) + 50;
T_an(xnode)

plot(xnode,T_ap,";4 Nodos;");
hold on;
grid on;

vec_x=linspace(L1,L2,cant_nodos*cant_nodos);
vec_y=T_an(vec_x);

plot(vec_x,vec_y,";Exacta(nodos^2);");
