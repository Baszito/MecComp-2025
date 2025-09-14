#PRUEBA DIF FIN
#ESTACIONARIO
clc;
clear all;
#xnode=[0 1/3 2/3 1];
#model=[2 0 0 0 100];
#cb=[1 10 0;1 50 0];
#et=0;

#[T_dir]=difFinitas(xnode,model,cb,et)


#xnode=linspace(0,2,6)
#model=[1 1 0 0 0];
#cb=[1 100 0;2 0 0];
#et=0;

#[T_neu]=difFinitas(xnode,model,cb,et)

#xnode=linspace(0,1,96)
#model.rho=0;
#model.cp=0;
#model.c=1;
#model.k=1;
#model.Q=50*ones(length(xnode),1);
#cb=[1 10 0;3 0.2 50];
#et=0;

#[T_rob]=difFinitas(xnode,model,cb,et)

#xnode=linspace(1,5,256);
#model.rho=0;
#model.cp=0;
#model.c=0;
#model.k=1;
#model.Q=100.*(xnode-3).^2;
#model.Q
#cb=[2 2 0;1 0 0];
#et=0;

#[T_neu]=difFinitas(xnode,model,cb,et);

xnode=linspace(5,10,16);
model.rho=1;
model.cp=1;
model.c=0;
model.k=2;
model.Q=(xnode).^3;
cb=[3 2 100;1 50 0];
et=1;

[T_rob]=difFinitas(xnode,model,cb,et);
