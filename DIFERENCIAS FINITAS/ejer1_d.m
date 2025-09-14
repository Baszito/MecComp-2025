#ejercicio 1 - D)
#---------DATOS--------
x=linspace(0,1,96);
n=length(x);


#---------Diferencias centradas con nodo ficticio--------#
T_an=@(x) (-36.6897.*exp(-x) - 3.3103.*exp(x) + 50);


#-------GRAFICO
figure;
vec_y = T_an(x);


plot(x,vec_y,';T analitico;');
