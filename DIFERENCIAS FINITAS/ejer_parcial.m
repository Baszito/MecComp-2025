clc;
clear all;
#Primero, expandi en series de taylor y obtene los coeficientes :
matriz=[1 1 1;-3 0 1/2;9/2 0 1/8];
#Segundo, arma los coef. ACORDATE DE PASAR EL QUE DIVIDE AL COEFICIENTE (tipo el 2!, 3! y asi) PARA EL OTRO LADO
b=[0 1 0];
b=b';
matriz\b


