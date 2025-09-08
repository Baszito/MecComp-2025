clc; clear all;

#Inicio y Final
L1 = 0;
L2 = 2;

#Constantes

pCp = 0;
k=1;
c=1;
G=0;

#Numero
N = 10;
deltax=(L2-L1)/(N-1);

#Condiciones de borde
T1=100;
q=0;

# A*phi_i+1 + B*phi_i + C*phi_i-1 = 2*deltax*G = D

A = pCp-(2*k/deltax);

B = (4*k/deltax) + (2*c*deltax);

C = -pCp - (2*k/deltax);

D = 2*deltax*G;


m = N - 2;
matriz = spdiags([A*ones(m,1), B*ones(m,1), C*ones(m,1)], [0, 1, 2], m, N);
matriz = full(matriz); % Convertir a matriz densa


matriz=[[1,zeros(1,N-1)];matriz];
ultima_fila = [zeros(1,N-2),[A+C,B]];
matriz = [matriz;ultima_fila];

x = matriz \ [T1,ones(1,N-2)*D,2*deltax*G-2*A*deltax*q/k]';

sol_analitica = @(x) ( 100*(e.^(-1*x)).*( (e.^(2*x)) + (e^4) ) ) / ( 1 + (e^4) );

vec_x = linspace(L1,L2,N);
vec_y = sol_analitica(vec_x)

plot(vec_x,x);
hold on;
grid on;
plot(vec_x,vec_y);
