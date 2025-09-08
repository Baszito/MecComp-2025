#ejercicio 1 - F)
clc;
clear;
#---------DATOS--------
k=2;
H=2;
c=2;
dx=1/3;#Cant de puntos que se quiere
L1=0;
L2=1;
pcp=2;
x=L1:dx:L2;
n=length(x);
Q=75*ones(length(x));
Qinf=10;
K=zeros(n);
b=zeros(n,1);

#TEMPORAL
T1=0;
T2=10;
dt=0.25*dx*dx/k;
Ttotal=T1:dt:T2;


#ARMADO MATRIZ K
#FUNCIONALIDAD FIJA
#NO VA SIGNO MENOS,OJO
fila = [k/dx^2 -(2*k/dx^2+c) k/dx^2];
  for i=2:(n-1)
    K(i,i-1:i+1)=fila;
    b(i)=Q(i);
end
K(1,:)=zeros(n,1);
K(n,n)=(-2*k/dx^2 -2*H/dx - c);
K(n,n-1)=2*k/dx^2;
K(1,1)=1;

b(n)=Q(end)+(2*H/dx)*Qinf;
b(1)=0;

#condiciones iniciales nulas
#T(nodo,tiempo)
IDENTIDAD=eye(n);
T_ap=zeros(length(x),length(Ttotal));
T_ap(1,1)=0;#CONDICION DIRICHLET

T_an=((-5/4).*e.^-(x+1)).*((e.^x)-1).*(11.*e.^x+(11-30*e));
ejes=[L1 L2 -0 20];

for i=1:length(Ttotal)
  figure(1)
  plot(x,T_ap(:,i),'*-')
  title(sprintf('t = %5.3f',i*dt));
  grid on;
  grid minor;
  axis(ejes);
  pause(0.01);
  #FORWARD
  #T_ap(:,i+1)=((dt/pcp)*(K*T_ap(:,i)+b))+T_ap(:,i);
  #T_ap(end,i+1)=50;#CONDICION DIRICHLET

  #BACKWARD
  #T_ap(:,i+1)=((pcp/dt)*IDENTIDAD-K)\(b+(pcp/dt)*T_ap(:,i));
  #T_ap(end,i+1)=50;#CONDICION DIRICHLET

  #CRANK-NICHOLSON
  T_ap(:,i+1)=((pcp/dt)*IDENTIDAD - 0.5*K)\(((0.5*K + (pcp/dt)*IDENTIDAD)*T_ap(:,i)) +b);
  T_ap(1,i+1)=0;#CONDICION DIRICHLET
endfor

hold on;
plot(x,T_an,'r');


