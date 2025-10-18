#---------------GUIA 1D - FEM---------------#
#Ejercicio 1 - C
#Aca No puedo despejar G, porque no es constante
#EJERCICIO EN COORDENADAS NATURALES
clc;
clear all;

#Datos:
pcp=0;
k=1;
c=0;
G=@(x) 100.*(x-3).^2;
H=0;

#Discretizacion
L1=1;
L2=5;
cant_nodos=4;
h=(L2-L1)/(cant_nodos-1);
xnode=L1:h:L2;
G=G(xnode);

#Condiciones de borde :
q_L1=2;
T_L2=0;

#Armado de matriz Ke

#Aporte Difusivo :
De=zeros(2,2);
De(1,1)=k/h;
De(1,2)=-k/h;
De(2,1)=-k/h;
De(2,2)=k/h;

#Aporte Reactivo #ESTA TA MAL
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


endfor

#imponemos condiciones de borde



#Armado vector F
#Resolviendo por cuadratura de gauss-legendre
F(1)+=167.901;
F(2)+=88.888;

F(2)+=9.8765;
F(3)+=9.8765;

F(3)+=88.888;
F(4)+=167.9012;
#Condicion Neumann :
F(1)+=-q_L1;

#Aporte Robin
K(cant_nodos,cant_nodos)=K(cant_nodos,cant_nodos)+H;

#Condicion Dirichlet :
K(cant_nodos,:)=zeros(1,cant_nodos);
K(cant_nodos,cant_nodos)=1;
F(cant_nodos)=T_L2;
F
#obtenemos la solucion aproximada
T_ap=K\F

#Solucion analitica, para verificar
T_an=@(y) (-25.*y.^4 +300.*y.^3 -1350.*y.^2 +1906.*y + 2345)/3;
T_an(xnode)

plot(xnode,T_ap,";4 Nodos;");
hold on;
grid on;

vec_x=linspace(L1,L2,cant_nodos*16);
vec_y=T_an(vec_x);

plot(vec_x,vec_y,";Exacta(4*h nodos);");
