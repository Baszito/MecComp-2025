#ejercicio 1 - b)
#---------DATOS--------

T_an=@(y) (-25.*y.^4 +300.*y.^3 -1350.*y.^2 +1906.*y + 2345)/3;

#-------GRAFICO
vec_x = linspace(1,5,256);
vec_y = T_an(vec_x);

