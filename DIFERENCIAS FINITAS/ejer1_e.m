#ejercicio 1 - E)
clc;
clear;
#---------DATOS--------
k=2;
H=2;
c=0;
dx=1/3;#Cant de puntos que se quiere
L1=5;
L2=10;
pcp=1;
x=L1:dx:L2;
n=length(x);
Q=@(x) x.^3;
Q=Q(x);
Qinf=100;
K=zeros(n);
b=zeros(n,1);

#TEMPORAL
T1=0;
T2=15;
dt=0.25*dx*dx/k;
Ttotal=T1:dt:T2;


#ARMADO MATRIZ K
#FUNCIONALIDAD FIJA
fila = [k/dx^2 -(2*k/dx^2+c) k/dx^2];
  for i=2:(n-1)
    K(i,i-1:i+1)=fila;
    b(i)=Q(i);
end
K(n,:)=zeros(n,1);
K(1,1)=(-2*k/dx^2 -2*H/dx - c);
K(1,2)=2*k/dx^2;
K(n,n)=1;

b(1)=Q(1)+(2*H/dx)*Qinf;
b(n)=50;

#condiciones iniciales nulas
#T(nodo,tiempo)
IDENTIDAD=eye(n);
T_ap=zeros(length(x),length(Ttotal));
T_ap(end,1)=50;#CONDICION DIRICHLET

T_an=(-x.^5)/40 + 1225.*x./3 -4600/3;
ejes=[L1 L2 -0 1000];

for i=1:length(Ttotal)
  figure(1)
  plot(x,T_ap(:,i),'*-')
  title(sprintf('t = %5.3f',i*dt));
  grid on;
  grid minor;
  axis(ejes);
  pause(0.01);
  #FORWARD
  T_ap(:,i+1)=((dt/pcp)*(K*T_ap(:,i)+b))+T_ap(:,i);
  T_ap(end,i+1)=50;#CONDICION DIRICHLET

  #BACKWARD
  #T_ap(:,i+1)=((pcp/dt)*IDENTIDAD-K)\(b+(pcp/dt)*T_ap(:,i));
  #T_ap(end,i+1)=50;#CONDICION DIRICHLET

  #CRANK-NICHOLSON
  #T_ap(:,i+1)=((pcp/dt)*IDENTIDAD - 0.5*K)\(((0.5*K + (pcp/dt)*IDENTIDAD)*T_ap(:,i)) +b);
  #T_ap(end,i+1)=50;#CONDICION DIRICHLET
endfor

hold on;
plot(x,T_an,'r');


