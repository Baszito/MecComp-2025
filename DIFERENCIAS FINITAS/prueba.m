#PRUEBA DIF FIN
#ESTACIONARIO
clc;
clear all;
xnode=[0 1/3 2/3 1];
model=[2 0 0 0 100];
cb=[1 10 0;1 50 0];
et=0;

[T_dir]=difFinitas(xnode,model,cb,et)


xnode=linspace(0,2,6)
model=[1 1 0 0 0];
cb=[1 100 0;2 0 0];
et=0;

[T_neu]=difFinitas(xnode,model,cb,et)

xnode=linspace(0,1,4)
model=[1 1 0 0 50];
cb=[1 10 0;3 0.2 50];
et=0;

[T_rob]=difFinitas(xnode,model,cb,et)
