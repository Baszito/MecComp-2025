clc;
clear all;
#Primero, expandi en series de taylor y obtene los coeficientes :
matriz=[1 1 1 1 1;0 1 -1 2 -2;0 1 1 4 4;0 1 -1 8 -8;0 1 1 16 16];
#Segundo, arma los coef. ACORDATE DE PASAR EL QUE DIVIDE AL COEFICIENTE (tipo el 2!, 3! y asi) PARA EL OTRO LADO
b=[0 0 0 6 0];
b=b';
matriz\b


